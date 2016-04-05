! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! Banner for checking maximum lengths (i.e. 32)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
! BEGIN DART PREPROCESS KIND LIST
! AMSRE_BRIGHTNESS_T,             KIND_BRIGHTNESS_TEMPERATURE
! END DART PREPROCESS KIND LIST

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_brightnessT_mod, only : get_amsre_metadata, &
!                                      read_amsre_metadata, &
!                                      write_amsre_metadata, &
!                                      interactive_amsre_metadata
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(AMSRE_BRIGHTNESS_T)
!     ! This operator reads information from a CLM history file, whose contents
!     ! are not modified by the assimilation process. There is no point doing
!     ! the posterior call as it returns the same value as the prior.
!     ! Need to pass metadata and ensemble index to interpolate. Terrible.
!     if (isprior) then
!        call interpolate(state, location, KIND_BRIGHTNESS_TEMPERATURE, obs_val, istatus, &
!        (/ get_amsre_metadata(obs_def%key), real(ens_index,r8) /) )
!     else
!        obs_val = MISSING_R8
!        istatus = 1
!     endif
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!  case(AMSRE_BRIGHTNESS_T)
!     call read_amsre_metadata(obs_def%key, key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!  case(AMSRE_BRIGHTNESS_T)
!     call write_amsre_metadata(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!  case(AMSRE_BRIGHTNESS_T)
!     call interactive_amsre_metadata(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_brightnessT_mod

! Prototype code for brightness temperature observation support.
! These observations have some extra metadata.

use        types_mod, only : r4, r8, MISSING_R8, MISSING_I
                            
use    utilities_mod, only : register_module, E_ERR, E_MSG, &
                             error_handler, ascii_file_format

use typesizes
use netcdf

implicit none
private

public ::  set_amsre_metadata, &
           get_amsre_metadata, &
          read_amsre_metadata, &
         write_amsre_metadata, &
   interactive_amsre_metadata

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical            :: module_initialized = .false.
character(len=512) :: string1, string2

!----------------------------------------------------------------------
! Metadata for AMSR-E observations.
!----------------------------------------------------------------------

type amsre_metadata
   private
   real(r8)  :: frequency    ! gHz
   real(r8)  :: footprint    ! km^2  'support' of the observation
   character :: polarization ! vertical or horizontal
   integer   :: landcovercode
end type amsre_metadata

type(amsre_metadata), allocatable, dimension(:) :: observation_metadata
type(amsre_metadata) :: missing_metadata

character(len=8), parameter :: AMSRESTRING = 'amsr-e'
character(len=8), parameter ::  IGBPSTRING = 'igbp'
real(r8),         parameter :: AMSRE_inc_angle = 55.0_r8 ! incidence angle (degrees)

integer, SAVE :: MAXamsrekey = 625000  ! 1 day, NH 25km, 1 freq, 1 pol, both passes
integer, SAVE ::    amsrekey = 0       ! useful length of metadata arrays

!----------------------------------------------------------------------
! no namelist items
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! This function name will be a problem if cosmos and amsrE are needed
! at the same time. refactor into single function in utilities_mod?
!----------------------------------------------------------------------

!interface interactive
!   module procedure interactive_real
!   module procedure interactive_int
!end interface


contains


!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of public routines for obs_def_brightnessT_mod
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine set_amsre_metadata(key, frequency, footprint, polarization, landcovercode)

! Fill the module storage metadata for a particular observation.

integer,   intent(out) :: key
real(r8),  intent(in)  :: frequency, footprint
character, intent(in)  :: polarization
integer,   intent(in)  :: landcovercode

if ( .not. module_initialized ) call initialize_module

amsrekey = amsrekey + 1  ! increase module storage used counter

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(amsrekey,'set_amsre_metadata')

key = amsrekey ! now that we know its legal

observation_metadata(key)%frequency     = frequency
observation_metadata(key)%footprint     = footprint
observation_metadata(key)%polarization  = polarization
observation_metadata(key)%landcovercode = landcovercode

end subroutine set_amsre_metadata


!======================================================================


function get_amsre_metadata(key)

! Query the metadata in module storage for a particular observation.
! encoding the polarization is required because Fortran is picky.
! realval = real('H',r8)   works fine, but the inverse does not:
! charval = char(int(realval)) ... FAILS ... so I decided to encode
! Horizontal polarizations as a positive value, and
! Vertical polarizations as a negative value.

integer,    intent(in) :: key
real(r8), dimension(5) :: get_amsre_metadata

if ( .not. module_initialized ) call initialize_module

! Make sure the desired key is within the length of the metadata arrays.
call key_within_range(key,'get_amsre_metadata')

get_amsre_metadata(1)  = real(observation_metadata(key)%landcovercode,r8)
get_amsre_metadata(2)  =      observation_metadata(key)%frequency
get_amsre_metadata(3)  =      observation_metadata(key)%footprint
if (observation_metadata(key)%polarization == 'H') then
   get_amsre_metadata(4)  = 1.0_r8
else
   get_amsre_metadata(4)  = -1.0_r8
endif
get_amsre_metadata(5)  =      AMSRE_inc_angle

end function get_amsre_metadata


!======================================================================


 subroutine read_amsre_metadata(key,       obsID, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_amsre_metadata(obs_def%key, key, ifile, fform)
!
! This routine reads the metadata for neutron intensity observations.
!
integer,          intent(out)          :: key    ! index into local metadata
integer,          intent(in)           :: obsID
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temp variables
integer           :: ierr
character(len=9)  :: header

real(r8)          :: frequency, footprint
character         :: polarization
integer           :: landcovercode

if ( .not. module_initialized ) call initialize_module

! write(string2,*)'observation #',obsID

if ( ascii_file_format(fform) ) then

   read(ifile, *, iostat=ierr) header
   call check_iostat(ierr,'read_amsre_metadata','header',string2)
   if (trim(header) /= trim(AMSRESTRING)) then
      write(string1,*)"Expected Tb header ["//AMSRESTRING//"] in input file, got ["//header//"]"
      call error_handler(E_ERR,'read_amsre_metadata',string1,source,revision,revdate,text2=string2)
   endif
   read(ifile, *, iostat=ierr) frequency, footprint, polarization
   call check_iostat(ierr,'read_amsre_metadata','freq foot polar',string2)

   read(ifile, *, iostat=ierr) header
   if (trim(header) /= trim(IGBPSTRING)) then
      write(string1,*)"Expected IGBP header ["//IGBPSTRING//"] in input file, got ["//header//"]"
      call error_handler(E_ERR,'read_amsre_metadata',string1,source,revision,revdate,text2=string2)
   endif
   read(ifile, *, iostat=ierr) landcovercode
   call check_iostat(ierr,'read_amsre_metadata','landcovercode',string2)

else

   read(ifile, iostat=ierr) header
   call check_iostat(ierr,'read_amsre_metadata','header',string2)
   if (trim(header) /= trim(AMSRESTRING)) then
      write(string1,*)"Expected Tb header ["//AMSRESTRING//"] in input file, got ["//header//"]"
      call error_handler(E_ERR,'read_amsre_metadata',string1,source,revision,revdate,text2=string2)
   endif
   read(ifile, iostat=ierr) frequency, footprint, polarization
   call check_iostat(ierr,'read_amsre_metadata','freq foot polar',string2)

   read(ifile, iostat=ierr) header
   if (trim(header) /= trim(IGBPSTRING)) then
      write(string1,*)"Expected IGBP header ["//IGBPSTRING//"] in input file, got ["//header//"]"
      call error_handler(E_ERR,'read_amsre_metadata',string1,source,revision,revdate,text2=string2)
   endif
   read(ifile, iostat=ierr) landcovercode
   call check_iostat(ierr,'read_amsre_metadata','landcovercode',string2)

endif

! Store the metadata in module storage and record the new length of the metadata arrays.
call set_amsre_metadata(key, frequency, footprint, polarization, landcovercode)

! The new 'key' is returned.

end subroutine read_amsre_metadata


!======================================================================


 subroutine write_amsre_metadata(key, ifile, fform)
!----------------------------------------------------------------------
! writes the metadata for AMSR-E brightness temperature observations.

integer,           intent(in)           :: key
integer,           intent(in)           :: ifile
character(len=*),  intent(in), optional :: fform


real(r8), dimension(5) :: metadata
real(r8)  :: frequency, footprint, incidence_angle
character :: polarization
integer   :: landcovercode


if ( .not. module_initialized ) call initialize_module

! given the index into the local metadata arrays - retrieve
! the metadata for this particular observation. All AMSR-E observations
! have the same incidence angle, so that is not written.

metadata = get_amsre_metadata(key)

landcovercode   =      int(metadata(1))
frequency       =          metadata(2)
footprint       =          metadata(3)
if (metadata(4) > 0.0_r8 ) then
   polarization = 'H'
else
   polarization = 'V'
endif
incidence_angle =          metadata(5)

if ( ascii_file_format(fform)) then
   write(ifile, *) trim(AMSRESTRING)
   write(ifile, '(f8.2,1x,f8.2,1x,A)') frequency, footprint, polarization
   write(ifile, *) trim(IGBPSTRING)
   write(ifile, '(i8)') landcovercode
else
   write(ifile) trim(AMSRESTRING)
   write(ifile) frequency, footprint, polarization
   write(ifile) trim(IGBPSTRING)
   write(ifile) landcovercode
endif

end subroutine write_amsre_metadata


!======================================================================
! The AMSR-E Level-2A product (AE_L2A) contains brightness temperatures at
! 6.9 GHz, 10.7 GHz, 18.7 GHz, 23.8 GHz, 36.5 GHz, and 89.0 GHz.
! Data are resampled to be spatially consistent, and therefore are available
! at a variety of resolutions that correspond to the footprint sizes of the
! observations such as 56 km, 38 km, 21 km, 12 km, and 5.4 km, respectively.
! Each swath is packaged with associated geolocation fields.
! Data are stored in Hierarchical Data Format - Earth Observing System (HDF-EOS)
! format and are available from 1 June 2002 to 4 October 2011 via FTP.

subroutine interactive_amsre_metadata( key )

integer, intent(out) :: key

real(r8)  :: frequency, footprint
character :: polarization
integer   :: landcovercode

if ( .not. module_initialized ) call initialize_module

! Prompt for input for the required metadata

frequency     = interactive_real('"frequency"    [GHz]',minvalue = 0.0_r8)
footprint     = interactive_real('"footprint"    [km]' ,minvalue = 0.0_r8)
polarization  = interactive_char('"polarization" [H,V]')
landcovercode = interactive_int( '"IGBP land cover code" [??]',minvalue = 0)

call set_amsre_metadata(key, frequency, footprint, polarization, landcovercode)

end subroutine interactive_amsre_metadata


!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of internal routines for obs_def_brightnessT_mod
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine initialize_module

! Called once to set values and allocate space.

integer :: iunit, io

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! Allocate the module array to store the key information.
missing_metadata%frequency     = MISSING_R8
missing_metadata%footprint     = MISSING_R8
missing_metadata%polarization  = 'x'
missing_metadata%landcovercode = MISSING_I

allocate( observation_metadata(MAXamsrekey) )
amsrekey = 0
observation_metadata(:) = missing_metadata

end subroutine initialize_module


!======================================================================


function interactive_real(str1,minvalue,maxvalue)
real(r8)                       :: interactive_real
character(len=*),   intent(in) :: str1
real(r8), optional, intent(in) :: minvalue
real(r8), optional, intent(in) :: maxvalue

integer :: i

interactive_real = MISSING_R8

! Prompt with a minimum amount of error checking

if     (present(minvalue) .and. present(maxvalue)) then

   interactive_real = minvalue - 1.0_r8
   MINMAXLOOP : do i = 1,10
      if ((interactive_real >= minvalue) .and. (interactive_real <= maxvalue)) exit MINMAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_real
   end do MINMAXLOOP

elseif (present(minvalue)) then

   interactive_real = minvalue - 1.0_r8
   MINLOOP : do i=1,10
      if (interactive_real >= minvalue) exit MINLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_real
   end do MINLOOP

elseif (present(maxvalue)) then

   interactive_real = maxvalue + 1.0_r8
   MAXLOOP : do i=1,10
      if (interactive_real <= maxvalue) exit MAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_real
   end do MAXLOOP

else ! anything goes ... cannot check
      write(*, *) 'Enter '//str1
      read( *, *) interactive_real
endif

end function interactive_real


!======================================================================


function interactive_char(str1)
character                      :: interactive_char
character(len=*),   intent(in) :: str1

interactive_char = ''

write(*, *) 'Enter '//str1
read( *, *) interactive_char

end function interactive_char


!======================================================================


function interactive_int(str1,minvalue,maxvalue)
integer                       :: interactive_int
character(len=*),  intent(in) :: str1
integer, optional, intent(in) :: minvalue
integer, optional, intent(in) :: maxvalue

integer :: i

interactive_int = MISSING_I

! Prompt with a minimum amount of error checking

if     (present(minvalue) .and. present(maxvalue)) then

   interactive_int = minvalue - 1
   MINMAXLOOP : do i = 1,10
      if ((interactive_int >= minvalue) .and. (interactive_int <= maxvalue)) exit MINMAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_int
   end do MINMAXLOOP

elseif (present(minvalue)) then

   interactive_int = minvalue - 1
   MINLOOP : do i=1,10
      if (interactive_int >= minvalue) exit MINLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_int
   end do MINLOOP

elseif (present(maxvalue)) then

   interactive_int = maxvalue + 1
   MAXLOOP : do i=1,10
      if (interactive_int <= maxvalue) exit MAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_int
   end do MAXLOOP

else ! anything goes ... cannot check
      write(*, *) 'Enter '//str1
      read( *, *) interactive_int
endif

end function interactive_int


!======================================================================


subroutine key_within_range(key, routine)

! Make sure we are addressing within the metadata arrays

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

! fine -- no problem.
if ((key > 0) .and. (key <= amsrekey)) return

! Bad news. Tell the user.
write(string1, *) 'key (',key,') not within known range ( 1,', amsrekey,')'
call error_handler(E_ERR,routine,string1,source,revision,revdate)

end subroutine key_within_range


!======================================================================


subroutine grow_metadata(key, routine)
!----------------------------------------------------------------------
! If the allocatable metadata arrays are not big enough ... try again

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

integer :: orglength
type(amsre_metadata), allocatable, dimension(:) :: safe_metadata

! fine -- no problem.
if ((key > 0) .and. (key <= MAXamsrekey)) return

orglength   =     MAXamsrekey
MAXamsrekey = 2 * orglength

! Check for some error conditions.
if (key < 1) then
   write(string1, *) 'key (',key,') must be >= 1'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
elseif (key >= 2*MAXamsrekey) then
   write(string1, *) 'key (',key,') really unexpected.'
   write(string2, *) 'doubling storage will not help.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate, &
                      text2=string2)
endif

! News. Tell the user we are increasing storage.
write(string1, *) 'key (',key,') exceeds Nmax_amsre_Tb (',orglength,')'
write(string2, *) 'Increasing Nmax_amsre_Tb to ',MAXamsrekey
call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)

allocate(safe_metadata(orglength))
safe_metadata(:) = observation_metadata(:)

deallocate(observation_metadata)
  allocate(observation_metadata(MAXamsrekey))

observation_metadata(1:orglength)             = safe_metadata(:)
observation_metadata(orglength+1:MAXamsrekey) = missing_metadata

deallocate(safe_metadata)

end subroutine grow_metadata


!======================================================================


subroutine check_iostat(istat, routine, varname, msgstring)

integer,          intent(in) :: istat
character(len=*), intent(in) :: routine
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: msgstring

if ( istat /= 0 ) then
   write(string1,*)'istat should be 0 but is ',istat,' for '//varname
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=msgstring)
end if

end subroutine check_iostat


!======================================================================


end module obs_def_brightnessT_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

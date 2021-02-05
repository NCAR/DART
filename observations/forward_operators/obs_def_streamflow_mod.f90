! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

!----------------------------------------------------------------------
! This module provides support for observations of streamflow.
!
!  OBS            1
!           -1           2          -1
! obdef
! loc3d
!      4.188790204786391        0.6986026390077337         137.7093918137252      3
! kind
!            1
! gageID_linkID
!    <a chararacter string ID for the gage>    <numeric link ID>
!      0          0
!    4.0000000000000000
!----------------------------------------------------------------------

! keep in mind that fortran allows only 31 characters in parameter
! definitions (which is what this string is going to be used for).
! if the platform name gets longer than 5 chars, consider going
! to something like xxx_TOTAL_PRECIP_WATER to give you room to
! put in more descriptive platform names.

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx-yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy-

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! STREAM_FLOW,                    QTY_STREAM_FLOW
! DEEP_GROUNDWATER_LEVEL,         QTY_DEEP_GROUNDWATER_LEVEL,    COMMON_CODE
! SKIN_TEMPERATURE,               QTY_SKIN_TEMPERATURE,          COMMON_CODE
! END DART PREPROCESS TYPE DEFINITIONS


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_streamflow_mod, only : read_streamflow_metadata, &
!                                     write_streamflow_metadata, &
!                               interactive_streamflow_metadata
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! The get_expected_obs_from_def for STREAM_FLOW is kinda dumb at the moment, but
! it expected that there will be a more complicated one soon.

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(STREAM_FLOW)
!     call interpolate(state_handle, ens_size, location, QTY_STREAM_FLOW, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(STREAM_FLOW)
!      call read_streamflow_metadata(obs_def%key, key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(STREAM_FLOW)
!      call write_streamflow_metadata(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(STREAM_FLOW)
!      call interactive_streamflow_metadata(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_streamflow_mod

use             types_mod, only : i4, r8, metadatalength

use         utilities_mod, only : register_module, error_handler, &
                                  E_ERR, E_MSG, &
                                  ascii_file_format

implicit none
private

public ::      set_streamflow_metadata, &
               get_streamflow_metadata, &
              read_streamflow_metadata, &
             write_streamflow_metadata, &
       interactive_streamflow_metadata, &
                       missing_gage_ID, &
                       missing_link_ID

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "obs_def_streamflow_mod.f90 $"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2
logical, save      :: module_initialized = .false.

! Metadata for streamflow observations.

integer, parameter :: IDLength = 15

type gage_metadata
   private
   character(len=IDLength) :: gage_ID
   integer(i4)             :: link_ID
end type gage_metadata

type(gage_metadata), allocatable, dimension(:) :: observation_metadata
type(gage_metadata) :: missing_metadata
character(len=*), parameter :: STREAMFLOWSTRING = 'gageID_linkID'

integer(i4) :: missing_link_ID = -1_i4
character(len=15) :: missing_gage_ID = 'UNKNOWN'

logical :: debug = .TRUE.
integer :: MAXstreamkey = 24*366  ! one year of hourly data - to start
integer ::    streamkey = 0       ! useful length of metadata arrays

!----------------------------------------------------------------------------
contains
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
!>

subroutine initialize_module

call register_module(source, revision, revdate)

module_initialized = .true.

missing_metadata%gage_ID = missing_gage_ID
missing_metadata%link_ID = missing_link_ID

allocate(observation_metadata(MAXstreamkey))

observation_metadata(:) = missing_metadata

end subroutine initialize_module


!----------------------------------------------------------------------
!> Fill the module storage metadata for a particular observation.

subroutine set_streamflow_metadata(key, gageID, linkID)

integer,          intent(out) :: key
character(len=*), intent(in)  :: gageID
integer(i4),      intent(in)  :: linkID

if ( .not. module_initialized ) call initialize_module

streamkey = streamkey + 1  ! increase module storage used counter

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(streamkey,'set_streamflow_metadata')

key = streamkey ! now that we know its legal

observation_metadata(key)%gage_ID = adjustl(gageID)
observation_metadata(key)%link_ID = linkID

end subroutine set_streamflow_metadata


!----------------------------------------------------------------------
!> Query the metadata in module storage for a particular observation.
!> This can be useful for post-processing routines, etc.

subroutine get_streamflow_metadata(key, gageID, linkID)

integer,          intent(in)  :: key
character(len=*), intent(out) :: gageID
integer(i4),      intent(out) :: linkID

if ( .not. module_initialized ) call initialize_module

! Make sure the desired key is within the length of the metadata arrays.
call key_within_range(key,'get_streamflow_metadata')

gageID = observation_metadata(key)%gage_ID
linkID = observation_metadata(key)%link_ID

end subroutine get_streamflow_metadata


!----------------------------------------------------------------------
!> This routine reads the metadata for streamflow observations.

subroutine read_streamflow_metadata(key, obsID, ifile, fform)
integer,          intent(out)          :: key    ! index into local metadata
integer,          intent(in)           :: obsID
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temp variables
logical           :: is_asciifile
integer           :: ierr
character(len=32) :: header
integer           :: oldkey
character(len=IDLength) :: gageID
integer(i4)             :: linkID

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

write(string2,*)'observation #',obsID

if ( is_asciifile ) then
   read(ifile, *, iostat=ierr) header
   call check_iostat(ierr,'read_streamflow_metadata','header',string2)
   if (header /= STREAMFLOWSTRING) then
       write(string1,*)"Expected streamflow header ["//STREAMFLOWSTRING//"] in input file, got ["//trim(header)//"]"
       call error_handler(E_ERR, 'read_streamflow_metadata', string1, source, revision, revdate, text2=string2)
   endif
   read(ifile, *, iostat=ierr) gageID, linkID
   call check_iostat(ierr,'read_streamflow_metadata','gageID_linkID',string2)
   read(ifile, *, iostat=ierr) oldkey
   call check_iostat(ierr,'read_streamflow_metadata','oldkey',string2)
else
   read(ifile, iostat=ierr) header
   call  check_iostat(ierr,'read_streamflow_metadata','header',string2)
   if (header /= STREAMFLOWSTRING) then
       write(string1,*)"Expected streamflow header ["//STREAMFLOWSTRING//"] in input file, got ["//trim(header)//"]"
       call error_handler(E_ERR, 'read_streamflow_metadata', string1, source, revision, revdate, text2=string2)
   endif
   read(ifile, iostat=ierr) gageID, linkID
   call check_iostat(ierr,'read_streamflow_metadata','gageID_linkID',string2)
   read(ifile, *, iostat=ierr) oldkey
   call check_iostat(ierr,'read_streamflow_metadata','oldkey',string2)
endif

! The oldkey is thrown away.

! Store the metadata in module storage and record the new length of the metadata arrays.
call set_streamflow_metadata(key, gageID, linkID)

! The new 'key' is returned.

end subroutine read_streamflow_metadata


!----------------------------------------------------------------------
!> writes the metadata for streamflow observations.

subroutine write_streamflow_metadata(key, ifile, fform)

integer,           intent(in)           :: key
integer,           intent(in)           :: ifile
character(len=*),  intent(in), optional :: fform

logical  :: is_asciifile
character(len=IDLength) :: gageID
integer(i4)             :: linkID

if ( .not. module_initialized ) call initialize_module

! given the index into the local metadata arrays - retrieve
! the metadata for this particular observation.

call get_streamflow_metadata(key, gageID, linkID)

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   write(ifile, *) trim(STREAMFLOWSTRING)
   write(ifile, *) gageID, linkID
   write(ifile, *) key
else
   write(ifile   ) trim(STREAMFLOWSTRING)
   write(ifile   ) gageID, linkID
   write(ifile   ) key
endif

end subroutine write_streamflow_metadata


!----------------------------------------------------------------------
!>

subroutine interactive_streamflow_metadata(key)

integer, intent(out) :: key

character(len=metadatalength) :: userinput
character(len=IDLength) :: gageID
integer(i4)             :: linkID
integer :: i, io

if ( .not. module_initialized ) call initialize_module

gageID = ''
linkID = -1

! Prompt for input for the required metadata

GAGELOOP : do i = 1,10
   write(string1,'("Enter up to ",i2," character")') IDLength
   write(*, *) trim(string1)//' "gage ID" string (NHD Gage Event ID ...)'
   read(*,*,iostat=io) userinput
   if (io == 0 .and. len_trim(userinput) <= IDLength) then
      gageID = trim(userinput)
      exit GAGELOOP
   endif
   if (io == 0 .and. len_trim(userinput)  > IDLength) &
      write(*,*)'too many characters (',len_trim(userinput),'), can only be ', IDLength
enddo GAGELOOP

LINKLOOP : do i = 1,10
   write(*, *) 'Enter numeric "Link ID" (NHDFlowline_network COMID)'
   read(*,*,iostat=io) linkID
   if (io == 0) exit LINKLOOP
enddo LINKLOOP

if (gageID /= '' .and. linkID > 0) then
   call set_streamflow_metadata(key, gageID, linkID)
else
   call error_handler(E_ERR,'interactive_streamflow_metadata', &
             'Unable to parse metadata.', source, revision, revdate)
endif

end subroutine interactive_streamflow_metadata


!----------------------------------------------------------------------
!>

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


!----------------------------------------------------------------------
!> Make sure we are addressing within the metadata arrays

subroutine key_within_range(key, routine)

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

! fine -- no problem.
if ((key > 0) .and. (key <= streamkey)) return

! Bad news. Tell the user.
write(string1, *) 'key (',key,') not within known range ( 1,', streamkey,')'
call error_handler(E_ERR,routine,string1,source,revision,revdate)

end subroutine key_within_range


!----------------------------------------------------------------------
!> If the allocatable metadata arrays are not big enough ... enlarge them

subroutine grow_metadata(key, routine)

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

integer :: orglength
type(gage_metadata), allocatable, dimension(:) :: safe_metadata

! fine -- no problem.
if ((key > 0) .and. (key <= MAXstreamkey)) return

orglength    =     MAXstreamkey
MAXstreamkey = 2 * orglength

! Check for some error conditions.
if (key < 1) then
   write(string1, *) 'key (',key,') must be >= 1'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
elseif (key >= 2*MAXstreamkey) then
   write(string1, *) 'key (',key,') really unexpected.'
   write(string2, *) 'doubling storage will not help.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate, &
                      text2=string2)
endif

! News. Tell the user we are increasing storage.
write(string1, *) 'key (',key,') exceeds MAXstreamkey (',orglength,')'
write(string2, *) 'Increasing MAXstreamkey to ',MAXstreamkey
call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)

allocate(safe_metadata(orglength))
safe_metadata(:) = observation_metadata(:)

deallocate(observation_metadata)
  allocate(observation_metadata(MAXstreamkey))

observation_metadata(1:orglength)              = safe_metadata(:)
observation_metadata(orglength+1:MAXstreamkey) = missing_metadata

deallocate(safe_metadata)

end subroutine grow_metadata



end module obs_def_streamflow_mod

! END DART PREPROCESS MODULE CODE


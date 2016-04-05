! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program TDT_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! TDT_to_obs - reads the TDT data as defined by no one.
! In a very forward-thinking move, the obs will have an associated IBGP type.
!
! original 10 March 2013   Tim Hoar   NCAR/IMAGe
!
! QC flags are defined as: 0 = BEST
!                          1 = OK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use          types_mod, only : r8, MISSING_R8, metadatalength

use      utilities_mod, only : initialize_utilities, finalize_utilities, &
                               register_module, error_handler, E_MSG, E_ERR, &
                               open_file, close_file, do_nml_file, do_nml_term, &
                               check_namelist_read, find_namelist_in_file, &
                               nmlfileunit, file_exist, nc_check, to_upper, &
                               find_textfile_dims

use   time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                               set_date, set_time, get_time, print_time, &
                               print_date, operator(-), operator(+), operator(>), &
                               operator(<), operator(==), operator(<=), operator(>=)

use       location_mod, only : location_type, set_location, VERTISHEIGHT

use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                               static_init_obs_sequence, init_obs, write_obs_seq, &
                               init_obs_sequence, get_num_obs, set_obs_def, &
                               set_copy_meta_data, set_qc_meta_data, &
                               set_obs_values, set_qc

use        obs_def_mod, only : obs_def_type, set_obs_def_kind, &
                               set_obs_def_location, set_obs_def_time, &
                               set_obs_def_error_variance, set_obs_def_key

use  obs_utilities_mod, only : add_obs_to_seq

! use    obs_def_TDT_mod, only : set_TDT_metadata

use       obs_kind_mod, only : SOIL_MOISTURE, SOIL_TEMPERATURE

use typesizes
use netcdf

implicit none

!-----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
!-----------------------------------------------------------------------

character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!-----------------------------------------------------------------------
! Namelist with default values
!-----------------------------------------------------------------------

character(len=128) :: text_input_file = 'corcounts.txt'
character(len=128) :: obs_out_file    = 'obs_seq.out'
integer            :: IBGB_class      = -1
real(r8)           :: latitude        = MISSING_R8
real(r8)           :: longitude       = MISSING_R8
real(r8)           :: maxgoodqc       = 99
logical            :: verbose         = .false.

namelist /TDT_to_obs_nml/ text_input_file, &
   obs_out_file, IBGB_class, latitude, longitude, maxgoodqc, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

character(len=256)      :: input_line, string1, string2, string3
integer                 :: iline
logical                 :: first_obs
integer                 :: oday, osec, rcio, iunit
integer                 :: num_copies, num_qc, max_obs
real(r8)                :: oerr, qc
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time
integer, parameter      :: LEVEL2QC = 1
integer                 :: obsindx


! The TDTdata type holds all the information from the csv file
! from a TDT instrument. At the moment, three standard depths
! are all that are supported.

type TDTdata
   integer               :: timeindex = 1
   character(len=20)     :: timestring = 'time'
   type(time_type)       :: time_obs
   integer               :: day    = 0
   integer               :: second = 0
   real(r8)              :: qc     = 1
   real(r8),dimension(3) :: swc   ! soil water content
   real(r8),dimension(3) :: temp  ! soil temperature
end type TDTdata

type(TDTdata) :: TDT ! we only need one of these.

! The site_metadata type holds all the site-specific information

type site_metadata
   type(location_type)   :: location
   real(r8)              :: latitude  = MISSING_R8
   real(r8)              :: longitude = MISSING_R8
   real(r8)              :: elevation = MISSING_R8
   real(r8)              :: IGBG ! code for site classification
   real(r8),dimension(3) :: depth     = (/ -0.05, 0.20, 0.50 /)
end type site_metadata

type(site_metadata) :: TDT_metadata

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('TDT_to_obs')

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "TDT_to_obs_nml", iunit)
read(iunit, nml = TDT_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "TDT_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=TDT_to_obs_nml)
if (do_nml_term()) write(     *     , nml=TDT_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)
prev_time = set_time(0, 0)

! Bulletproof these ... TJH FIXME
TDT_metadata%latitude  = latitude
TDT_metadata%longitude = longitude
TDT_metadata%elevation = elevation
TDT_metadata%IGBG      = IGBG_class

if (verbose) print *, 'TDT site located at lat, lon, elev  =', &
                       TDT_metadata%latitude, &
                       TDT_metadata%longitude, &
                       TDT_metadata%elevation

! We need to know the maximum number of observations in the input file.
! Each line has info for a single observation we want (TDT neutron counts).
! Each observation in this series will have a single
! observation value, its standard deviation, and a quality control flag.  
! Initialize two empty observations - one to track location
! in observation sequence - the other is for the new observation.

call find_textfile_dims(text_input_file, max_obs)

if (max_obs < 0) then
   write (string1,*) '<'//trim(text_input_file)//'> does not exist.'
   call error_handler(E_ERR,'main', string1, source, revision, revdate)
elseif (max_obs < 2) then
   write (string1,*) trim(text_input_file)//' has no observation values in it.'
   call error_handler(E_ERR,'main', string1, source, revision, revdate)
endif

iunit = open_file(text_input_file, 'formatted', 'read')
if (verbose) print *, 'opened input file ' // trim(text_input_file)

num_copies = 1
num_qc     = 1
first_obs  = .true.

call static_init_obs_sequence()
call init_obs(        obs,      num_copies, num_qc)
call init_obs(        prev_obs, num_copies, num_qc)
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(  obs_seq, 1, 'TDT QC')

! The first line describes all the fields ... column headers, if you will
call decode_header(iunit)

obsloop: do iline = 2,max_obs

   ! read in entire text line into a buffer
   read(iunit,'(A)',iostat=rcio) input_line
   if (rcio < 0) exit obsloop
   if (rcio > 0) then
      write (string1,'(''Cannot read (error '',i3,'') line '',i8,'' in '',A)') &
                    rcio, iline, trim(text_input_file)
      call error_handler(E_ERR,'main', string1, source, revision, revdate)
   endif

   ! parse the line into the TDT structure (including the observation time)
   call stringparse(input_line, iline)

   if (iline <= 2) then
      write(*,*)''
      write(*,*)'Check of the first observation: (column,string)'
      write(*,*)TDT%timeindex, TDT%timestring
      call print_date(TDT%time_obs, 'observation date is')
      call print_time(TDT%time_obs, 'observation time is')
   endif

   if (verbose) call print_date(TDT%time_obs, 'obs time is')

   call get_time(TDT%time_obs, osec, oday)

   ! make an obs derived type, and then add it to the sequence
   ! QC flags are pulled out of thin air.
   ! At some point, set_TDT_metadata() will employ a lookup table
   ! to provide the IGBG classification for a particular lat/lon
   ! Presently, it is specified in the namelist.

   DEPTH : DO idepth=1,length(TDT%swc)
      if (TDT%swc(idepth) /= MISSING_R8) then
         swc  = TDT%swc(idepth)
         oerr = min(0.1_r8, swc*0.1_r8) ! bald-faced guess
         qc   = 1                       ! bald-faced guess
      !  call set_TDT_metadata(obsindx, TDT_metadata%IGBG) 
         call create_3d_obs(TDT_metadata%latitude, TDT_metadata%longitude, &
              TDT_metadata%idepth, VERTISHEIGHT, swc, SOIL_MOISTURE, &
              oerr, oday, osec, qc, obsindx, obs)
         call add_obs_to_seq(obs_seq, obs, TDT%time_obs, prev_obs, prev_time, first_obs)
      endif
      if (TDT%temp(idepth) /= MISSING_R8) then
         temp  = TDT%temp(idepth)
         oerr = min(0.1_r8, abs(temp)*0.1_r8) ! bald-faced guess
         qc   = 1                             ! bald-faced guess
      !  call set_TDT_metadata(obsindx, TDT_metadata%IGBG) 
         call create_3d_obs(TDT_metadata%latitude, TDT_metadata%longitude, &
              TDT_metadata%idepth, VERTISHEIGHT, temp, SOIL_TEMPERATURE, &
              oerr, oday, osec, qc, obsindx, obs)
         call add_obs_to_seq(obs_seq, obs, TDT%time_obs, prev_obs, prev_time, first_obs)
      endif
   enddo DEPTH

enddo obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.  
!
!   NOTE: assumes the code is using the threed_sphere locations module, 
!         that the observation has a single data value and a single
!         qc value.
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    key   - index to metadata in obs_def_TDT_mod arrays
!    obs   - observation type
!
!    extended from the observations/utilities/obs_utilities_mod.f90 v 5601
!    to support the extra metadata -- TJH
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, key, obs)
integer,        intent(in)    :: okind, vkind, day, sec
real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
integer,        intent(in)    :: key
type(obs_type), intent(inout) :: obs

real(r8)              :: obs_val(1), qc_val(1)
type(obs_def_type)    :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def_key(obs_def, key)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

end subroutine create_3d_obs




subroutine decode_header(iunit)
! Reads the first line of the header and parses the information.
! I should break the line into words and match which word with each
! desired string. But not today ...
! FIXME ... decode the header ... do not assume ...
!
! head TDT_data-TIME-SWC5cm20cm50cm-TempC5cm20cm50cm.csv 
!        [      soil water content       ][     soil temperature     ]
!   time,       5cm,      20cm,      50cm,     5cm,     20cm,     50cm
! 734643,0.15780275,0.21564578,0.24085209,14.14314,12.974549,11.606505


integer, intent(in) :: iunit

read(iunit,'(A)',iostat=rcio) input_line
if (rcio /= 0) then
  write(string1,*)'Cannot parse header. Begins <',trim(input_line(1:40)),'>'
  call error_handler(E_ERR,'decode_header',string1, source, revision, revdate)
endif

call error_handler(E_MSG,'decode_header','hardcoding values for now ... dangerous', &
                     source, revision, revdate)

TDT%timeindex       = 1

end subroutine decode_header



subroutine stringparse(str1,linenum)
! The first two items are character strings, the rest can be read as reals.

character(len=*), intent(in) :: str1
integer         , intent(in) :: linenum

character(len=10) :: yyyymmdd
character(len=5)  :: tod
real(r8), dimension(9) :: values

values = MISSING_R8

left off here ....


read(str1,*,iostat=rcio) yyyymmdd,tod,values
if (rcio /= 0) then
  write(string1,*)'Cannot parse line',linenum,'. Begins <',trim(str1(1:30)),'>'
  call error_handler(E_ERR,'stringparse',string1, source, revision, revdate)
endif

read(yyyymmdd,'(i4,1x,i2,1x,i2)',iostat=rcio) TDT%year, TDT%month, TDT%day
if (rcio /= 0) then
  write(string1,*)'Cannot parse yyyymmdd on line ',linenum,' from <',yyyymmdd,'>'
  call error_handler(E_ERR,'stringparse',string1, source, revision, revdate)
endif

read(     tod,'(      i2,1x,i2)',iostat=rcio) TDT%hour, TDT%minute
if (rcio /= 0) then
  write(string1,*)'Cannot parse tod on line ',linenum,' from <',tod,'>'
  call error_handler(E_ERR,'stringparse',string1, source, revision, revdate)
endif

! Stuff what we want into the TDT structure

TDT%qc         = LEVEL2QC                           ! for all level 2 data
TDT%neutron    = values(TDT%neutronindex    - 2) ! minus the first two words
TDT%neutronstd = values(TDT%neutronstdindex - 2) ! minus the first two words
TDT%time_obs   = set_date(TDT%year, TDT%month, TDT%day, TDT%hour, TDT%minute, 0)

TDT%time_offset = 584755  ! 1601,1,1  in matlab-speak

end subroutine stringparse





end program TDT_to_obs

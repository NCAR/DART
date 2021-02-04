! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program prepbufr_to_obs

use  types_mod,        only : r8, rad2deg, PI
use     utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                              register_module, error_handler, &
                              timestamp, is_longitude_between, &
                              initialize_utilities, register_module,            &
                              do_output, nmlfileunit, do_nml_file, do_nml_term, &
                              finalize_utilities, E_ERR, E_MSG,  &
                              find_namelist_in_file, check_namelist_read
use  obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def, &
                              write_obs_def, destroy_obs_def, interactive_obs_def, &
                              copy_obs_def, set_obs_def_time, set_obs_def_type_of_obs, &
                              set_obs_def_error_variance, set_obs_def_location
use  time_manager_mod, only : time_type, operator(>), operator(<), operator(<=), &
                              operator(==), operator(/), operator(+), operator(-), &
                              set_date, set_calendar_type, get_time, get_date, &
                              set_time, GREGORIAN, print_date, print_time
use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs
use      location_mod, only : location_type, set_location, &
                              VERTISPRESSURE, VERTISSURFACE
use  obs_sequence_mod, only : init_obs_sequence, init_obs, obs_sequence_type, obs_type, &
                              set_copy_meta_data, set_qc_meta_data, write_obs_seq, &
                              static_init_obs_sequence, destroy_obs_sequence 
use      obs_kind_mod       ! all defined kinds
use obs_def_altimeter_mod, only: compute_altimeter


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(obs_sequence_type) :: seq

character(len=256) :: msgstring1, msgstring2
character(len=256) :: input_name, output_name
character(len=8  ) :: obsdate
integer :: iunit, io, ii, day1, kbeg, kend, obs_unit
integer :: fyear, fmonth, fday, fhour, fmin, fsec
integer :: hyear, hmonth, hday, hhour, hmin, hsec
integer :: gdays, gsecs, sdays, edays, dummy
logical :: daily, inc_midnight

type(time_type) :: start_time, end_time, window_width, one_sec
type(time_type) :: window_start, window_end, window_mid, window_half
type(time_type) :: hack, next


! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------
        
integer :: startyear = 2003, startmonth = 1, startday = 1, starthour = 0, startseconds = 0
integer ::   endyear = 2003,   endmonth = 1,   endday = 1,   endhour = 6,   endseconds = 0
integer :: windowdays = 0, windowhours = 6, windowseconds = 0
logical :: midnight_24 = .false.

integer :: max_num_obs = 2000000
character(len = 128) ::  input_filename_base    = 'temp_obs.'
character(len = 128) ::  input_filename_pattern = '(A,I4.4,3(I2.2))'  ! temp_obs.YYYYMMDDHH
character(len = 128) :: output_filename_base    = 'obs_seq.'
character(len = 128) :: output_filename_pattern = '(A,I4.4,3(I2.2))'  ! obs_seq.YYYYMMDDHH

logical :: select_obs = .false.

! FIXME:
! these aren't complete - there's a longer list in the prepbufr converter
logical :: ADPUPA = .false., AIRCAR = .false., AIRCFT = .false., &
           SATWND = .false., SATEMP = .false., SFCSHP = .false., &
           ADPSFC = .false.

! FIXME: check to see if these are complete
!  the prepbufr rec types are: 'P','Q','T','Z','U','V'
! what's a Z?  is it altimeter, height?
!------------------------------------------------------------------------------
logical :: obs_U  = .false., obs_V  = .false., obs_T  = .false. , &
           obs_PS = .false., obs_QV = .false.

real(r8) :: lon1 =   0.0_r8,  &   !  lower longitude bound
            lon2 = 360.0_r8,  &   !  upper longitude bound 
            lat1 = -90.0_r8,  &   !  lower latitude bound
            lat2 =  90.0_r8       !  upper latitude bound

logical  :: include_specific_humidity = .true.,  &
            include_relative_humidity = .false., &
            include_dewpoint          = .false.

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.
logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000
logical :: debug            = .false.

namelist /prepbufr_to_obs_nml/ &
        startyear, startmonth, startday, starthour, &
        endyear, endmonth, endday, endhour, &
        windowdays, windowhours, windowseconds, &
        input_filename_base, input_filename_pattern, &
        output_filename_base, output_filename_pattern, &
        max_num_obs, select_obs, midnight_24, &
        ADPUPA, AIRCAR, AIRCFT, SATEMP, SFCSHP, ADPSFC, SATWND, &
        obs_U, obs_V, obs_T, obs_PS, obs_QV, &
        lon1, lon2, lat1, lat2, &
        include_specific_humidity, include_relative_humidity, include_dewpoint, &
        print_timestamps, print_every_Nth, debug

! ----------------------------------------------------------------------
! Select observation types using NCEP categories (when select_obs is true).
!  ADPUPA: upper-air reports (mostly radiosonde plus few dropsonde, PIBAL)
!  AIRCFT: Conv. (AIREP, PIREP) and ASDAR aircraft reports
!  AIRCAR: ACARS sircraft reports
!  SATEMP: ATOVS retrived temperature
!  SFCSHP: SURFACE MARINE reports
!  ADPSFC: SURFACE LAND SYNOPTIC STATION reports
!  SATWND: Satellite derived wind reports
! ----------------------------------------------------------------------
! Select variables of U, V, T, QV, PS using the logicals:
!  obs_U   obs_V   obs_PS   obs_T   obs_QV  
! ----------------------------------------------------------------------

! start of executable program code


! Initialization and read namelist

call initialize_utilities('prepbufr_to_obs')
call register_module(source,revision,revdate)
call static_init_obs_sequence()

! Read the namelist entry
call find_namelist_in_file("input.nml", "prepbufr_to_obs_nml", iunit)
read(iunit, nml = prepbufr_to_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "prepbufr_to_obs_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=prepbufr_to_obs_nml)
if (do_nml_term()) write(     *     , nml=prepbufr_to_obs_nml)

! set observation calendar type and start/end dates
call set_calendar_type(GREGORIAN)

start_time = set_date(startyear, startmonth, startday, starthour, 0, startseconds)
  end_time = set_date(  endyear,   endmonth,   endday,   endhour, 0,   endseconds)

window_width = set_time(windowhours*3600 + windowseconds, windowdays)
if (window_width == set_time(0, 1)) then
   daily = .true.
else
   daily = .false.
endif
! FIXME: check to be sure width is even
window_half = window_width / 2
one_sec = set_time(1)

   call print_date(start_time,  ' Start time: ')
   call print_date(  end_time,  '   End time: ')
   call print_time(window_width,'Window size: ')   

! set up initial window before while loop below.
! the window start needs to be +1 second so no times
! are in more than a single window.
window_start = start_time 
window_end = window_start + window_width
window_mid = window_start + window_half
window_start = window_start + one_sec

! Loop through the time period of interest.
do while (window_start <= end_time) 
  
   ! output current window info
   call print_date(window_start, 'Window start time: ')
   call get_time(window_start, dummy, sdays)
   call print_date(window_end,   'Window   end time: ')
   call get_time(window_end, dummy, edays)
   
   ! FIXME: we could let user specify how to name the files
   ! with a window offset - 0 for start time, 1/2 for mid, 
   ! whole window width for end time...   because 6H files are
   ! typically named for the midpoint, but daily for 0Z and
   ! 3H files for the start time...   standards, so many.
   call print_date(window_mid, 'time for output file: ')
   
   call get_date(window_mid, fyear, fmonth, fday, fhour, fmin, fsec)
   call get_time(window_mid, gsecs, gdays)

   ! if this window crosses over a day boundary, add 24 to the
   ! end bin time in the construct_obs_sequence routine.  right now
   ! the ascii intermediate files only have hours and no dates in them.
   ! the current converter which creates them adds 24 to hours beyond
   ! midnight for files which are centered on midnight and include obs
   ! from 21Z one day to 3Z the next day.

   inc_midnight = (sdays /= edays)
   
   ! construct input and output filenames
   ! FIXME: make this consistent with prepbufr and the non-daily option:
   if (fhour == 0) then
      if (midnight_24) then
         hack = window_mid - set_time(0, 1)
         call get_date(hack, hyear, hmonth, hday, hhour, hmin, hsec)
         write( input_name,  input_filename_pattern) trim(input_filename_base), hyear, hmonth, hday, 24
      else
         hack = window_mid
         call get_date(hack, hyear, hmonth, hday, hhour, hmin, hsec)
         write( input_name,  input_filename_pattern) trim(input_filename_base), hyear, hmonth, hday, 0
      endif
   else
      write( input_name,  input_filename_pattern) trim(input_filename_base), fyear, fmonth, fday, fhour
   endif
   write(output_name, output_filename_pattern) trim(output_filename_base), fyear, fmonth, fday, fhour
   
   write(*,*) ' input filename: "'//trim( input_name)//'"'
   write(*,*) 'output filename: "'//trim(output_name)//'"'

   ! Initialize an obs_sequence 
      
   call init_obs_sequence(seq, 1, 1, max_num_obs)
      
   ! set meta data of obs_seq
      
   call set_copy_meta_data(seq, 1, 'NCEP BUFR observation')
   call set_qc_meta_data(seq, 1, 'NCEP QC index')

   write(*,*) ' opening main input file: "'//trim(input_name)//'"'
   obs_unit = open_file(input_name, form='formatted', action='read')
   
   call construct_obs_sequence(seq, obs_unit, gdays, window_start, window_end, inc_midnight)
   call close_file(obs_unit)

   ! construct input and output filenames
   ! FIXME: make this consistent with prepbufr and the non-daily option:
   next = window_mid + window_width
   call get_date(next, fyear, fmonth, fday, fhour, fmin, fsec)
   if (fhour == 0) then
      hack = next - set_time(0, 1)
      call get_date(hack, hyear, hmonth, hday, hhour, hmin, hsec)
      write( input_name,  input_filename_pattern) trim(input_filename_base), hyear, hmonth, hday, 24
   else
      write( input_name,  input_filename_pattern) trim(input_filename_base), fyear, fmonth, fday, fhour
   endif

   write(*,*) ' opening aux  input file: "'//trim(input_name)//'"'
   obs_unit = open_file(input_name, form='formatted', action='read')
   
   ! read the next available ascii intermediate file to collect any obs which are
   ! exactly equal to the ending timestamp.  this works for windows which are an even
   ! multiple of 6H (the times in the original prepbufr files).
   call construct_obs_sequence(seq, obs_unit, gdays, window_end, window_end, .false.)
   call close_file(obs_unit)
   
   ! output the sequence to a file
   call write_obs_seq(seq, output_name)
   
   ! release the memory of the seq.
   call destroy_obs_sequence(seq)
   

   window_start = window_end
   window_end = window_start + window_width
   window_mid = window_start + window_half
   window_start = window_start + one_sec

enddo

call error_handler(E_MSG,'prepbufr_to_obs','Finished successfully.',source,revision,revdate)
call finalize_utilities()

contains


!-------------------------------------------------


subroutine construct_obs_sequence (seq, iunit, gday, wstart, wend, wrap)
!------------------------------------------------------------------------------
!  this function is to convert NCEP decoded BUFR data to DART sequence format
!
!  FIXME: the prepbufr converter HAS to start putting a date as well as a time
!  in the ascii intermediate files, so i can jettison both the wrap and keep3
!  flags, and so we can have a single converter source (right now there is
!  prepbufr.f and prepbufr_03Z.f, which differ by about 3 lines.)
!

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: iunit
integer,                 intent(in)    :: gday
type(time_type),         intent(in)    :: wstart, wend
logical,                 intent(in)    :: wrap

type(obs_type) :: obs, prev_obs
integer :: i, num_copies, num_qc
integer :: seconds
integer :: hour, imin, sec
integer :: obs_num
integer :: syear, smonth, sday, shour, sminute, sseconds
integer :: eyear, emonth, eday, ehour, eminute, eseconds
type(time_type) :: current_day, time_obs, prev_time
real(r8) :: bin_beg, bin_end

! these can be initialized because they don't change
! for multiple calls.
integer, parameter :: num_fail_kinds = 8
integer :: iskip(num_fail_kinds)
integer, parameter :: fail_3Z        = 1
integer, parameter :: fail_timerange = 2
integer, parameter :: fail_badloc    = 3
integer, parameter :: fail_areabox   = 4
integer, parameter :: fail_badkind   = 5
integer, parameter :: fail_notwanted = 6
integer, parameter :: fail_badvert   = 7
integer, parameter :: fail_moisttype = 8
character(len=32) :: skip_reasons(num_fail_kinds) = (/ &
                     'time too early (exactly 03Z)    ', &
                     'time outside bin range          ', &
                     'bad observation location        ', &
                     'observation outside lat/lon box ', &
                     'unrecognized observation kind   ', &
                     'obs type not on select list     ', &
                     'bad vertical coordinate         ', &
                     'unwanted moisture type          ' /)


integer :: rday, read_counter, io
integer :: obs_prof, obs_kind, obs_kind_gen, which_vert, iqc, obstype, pc
real (r8) :: obs_err, lon, lat, lev, zob, time, rcount, zob2
real (r8) :: vloc, obs_value


character(len = 8 ) :: obsdate
character(len = 80) :: obsfile, label
character(len = 6 ) :: subset
logical :: pass, first_obs

call get_date(wstart, syear, smonth, sday, shour, sminute, sseconds)
call get_date(wend,   eyear, emonth, eday, ehour, eminute, eseconds)
if (debug) call print_date(wstart, 'window start')
if (debug) call print_date(wend,   'window   end')

bin_beg = shour + sminute / 60.0_r8 + sseconds / 3600.0_r8
bin_end = ehour + eminute / 60.0_r8 + eseconds / 3600.0_r8
if (bin_end < bin_beg) then
   bin_end = bin_end + 24.0_r8
endif
if (debug) print *, 'bin start, end: ', bin_beg, bin_end


! Initialize the obs variable

num_copies  = 1
num_qc      = 1

call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

obs_num = 0
iskip(:) = 0
first_obs = .true.

!  loop over all observations within the file
!------------------------------------------------------------------------------

read_counter = 1

obsloop:  do

   read(iunit,880,end=200,iostat=io) obs_err, lon, lat, lev, zob, zob2, rcount, time, &
                                        obstype, iqc, subset, pc
   if (io /= 0) then
      write(msgstring1,*)'read error was ',io,' for line ',read_counter
      call error_handler(E_ERR,'construct_obs_sequence', msgstring1, source, revision, revdate)
   endif

 880 format(f5.2,2f9.4,e12.5,f7.2,f7.2,f9.0,f7.3,i4,i2,1x,a6,i2)

   read_counter = read_counter + 1

   !write(*, 880) obs_err, lon, lat, lev, zob, zob2, rcount, time, &
   !                           obstype, iqc, subset, pc

!------------------------------------------------------------------------------
   !   A 'day' is from 03:01Z of one day through 03Z of the next.
   !   skip the observations at exact 03Z of the beginning of the day
   !   (obs at 03Z the next day have a time of 27.)
   ! (this should be deprecated eventually when we start putting dates
   ! as well as times in the intermediate ascii files.  for now, keep it.)
   if(time == 3.0_r8 .and. daily) then
      if (debug) write(*,*) 'invalid time.  hours = ', time
      iskip(fail_3Z) = iskip(fail_3Z) + 1
      cycle obsloop 
   endif 

   ! reject obs outside the current time window
   if(time < bin_beg .or. time > bin_end) then
      if (debug) write(*,*) 'invalid time.  hours = ', time
      iskip(fail_timerange) = iskip(fail_timerange) + 1
      cycle obsloop
   endif

   ! reject locations outside the valid values (shouldn't happen)
   if((lon > 360.0_r8) .or. (lon <   0.0_r8) .or.  &
      (lat >  90.0_r8) .or. (lat < -90.0_r8)) then
      if (debug) write(*,*) 'invalid location.  lon,lat = ', lon, lat
      iskip(fail_badloc) = iskip(fail_badloc) + 1
      cycle obsloop
   endif

   ! reject observations outside the bounding box
   if(lat < lat1 .or. lat > lat2 .or. & 
     .not. is_longitude_between(lon, lon1, lon2)) then
      if (debug) write(*,*) 'invalid location.  lon,lat = ', lon, lat
      iskip(fail_areabox) = iskip(fail_areabox) + 1
      cycle obsloop
   endif

   ! it seems when a machine writes out a number and then reads it back in
   ! you should at least get consistent round-off errors.  perhaps it is more
   ! complex than it seems, but on the ibm xlf compiler, values of exactly
   ! 90.0 and -90.0 are written out as (binary) doubles and when they are
   ! read back in and compared to 90 and -90, they are not equal anymore.
   ! the difference is in a very small but apparently significant digit,
   ! so we error out in our 3d sphere locations module.  this hack is a try
   ! at adjusting the numbers by a small enough value that the location is
   ! still the pole but isn't going to round out of range.
   if      (lat >=  89.9999_r8) then
     lat = lat - 1.0e-12_r8
     if (debug) write(*,*) 'lat adjusted down, now ', lat
   else if (lat <= -89.9999_r8) then
     lat = lat + 1.0e-12_r8
     if (debug) write(*,*) 'lat adjusted   up, now ', lat
   endif

   obs_prof = rcount/1000000


!   assign each observation the correct observation type
!------------------------------------------------------------------------------

   ! make sure we do not fall through the code below without setting
   ! a valid obs kind (e.g. the obstype is one not listed)
   obs_kind = -1

   if(obs_prof == 1) then
     obs_kind_gen = QTY_TEMPERATURE
     if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_TEMPERATURE
     if(obstype == 130 .or. obstype == 131) obs_kind = AIRCRAFT_TEMPERATURE
     if(obstype == 133                    ) obs_kind = ACARS_TEMPERATURE
     if(obstype == 161 .or. obstype == 163) obs_kind = ATOV_TEMPERATURE
     if(obstype == 171 .or. obstype == 173) obs_kind = ATOV_TEMPERATURE
     if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_TEMPERATURE
     if(obstype == 181 .or. obstype == 183) obs_kind = LAND_SFC_TEMPERATURE
   endif

   if(obs_prof == 5) then
     if ( zob2 == 0.0_r8 .and. include_specific_humidity ) then
       obs_kind_gen = QTY_SPECIFIC_HUMIDITY
       if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_SPECIFIC_HUMIDITY
       if(obstype == 130 .or. obstype == 131) obs_kind = AIRCRAFT_SPECIFIC_HUMIDITY
       if(obstype == 133                    ) obs_kind = ACARS_SPECIFIC_HUMIDITY
       if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_SPECIFIC_HUMIDITY
       if(obstype == 181 .or. obstype == 183) obs_kind = LAND_SFC_SPECIFIC_HUMIDITY
     else if ( zob2 == 1.0_r8 .and. include_relative_humidity ) then
       obs_kind_gen = QTY_RELATIVE_HUMIDITY
       if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_RELATIVE_HUMIDITY
       if(obstype == 130 .or. obstype == 131) obs_kind = AIRCRAFT_RELATIVE_HUMIDITY
       if(obstype == 133                    ) obs_kind = ACARS_RELATIVE_HUMIDITY
       if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_RELATIVE_HUMIDITY
       if(obstype == 181 .or. obstype == 183) obs_kind = LAND_SFC_RELATIVE_HUMIDITY
     else if ( zob2 == 2.0_r8 .and. include_dewpoint ) then
       obs_kind_gen = QTY_DEWPOINT
       if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_DEWPOINT
       if(obstype == 130 .or. obstype == 131) obs_kind = AIRCRAFT_DEWPOINT
       if(obstype == 133                    ) obs_kind = ACARS_DEWPOINT
       if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_DEWPOINT
       if(obstype == 181 .or. obstype == 183) obs_kind = LAND_SFC_DEWPOINT
     endif
   endif

   if(obs_prof == 3) then
     obs_kind_gen = QTY_SURFACE_PRESSURE
     if(obstype == 120 .or. obstype == 132) obs_kind = RADIOSONDE_SURFACE_ALTIMETER
     if(obstype == 180 .or. obstype == 182) obs_kind = MARINE_SFC_ALTIMETER
     if(obstype == 181                    ) obs_kind = LAND_SFC_ALTIMETER
   endif

   if(obs_prof == 2) then
     obs_kind_gen = QTY_U_WIND_COMPONENT
     if(obstype == 220 .or. obstype == 232) obs_kind = RADIOSONDE_U_WIND_COMPONENT
     if(obstype == 221                    ) obs_kind = RADIOSONDE_U_WIND_COMPONENT
     if(obstype == 230 .or. obstype == 231) obs_kind = AIRCRAFT_U_WIND_COMPONENT
     if(obstype == 233                    ) obs_kind = ACARS_U_WIND_COMPONENT
     if(obstype == 242 .or. obstype == 243) obs_kind = SAT_U_WIND_COMPONENT
     if(obstype == 245 .or. obstype == 246) obs_kind = SAT_U_WIND_COMPONENT
     if(obstype == 252 .or. obstype == 253) obs_kind = SAT_U_WIND_COMPONENT
     if(obstype == 255                    ) obs_kind = SAT_U_WIND_COMPONENT
     if(obstype == 280 .or. obstype == 282) obs_kind = MARINE_SFC_U_WIND_COMPONENT
     if(obstype == 281 .or. obstype == 284) obs_kind = LAND_SFC_U_WIND_COMPONENT
   endif

   if(obs_prof == 9) then
     obs_kind_gen = QTY_V_WIND_COMPONENT
     if(obstype == 220 .or. obstype == 232) obs_kind = RADIOSONDE_V_WIND_COMPONENT
     if(obstype == 221                    ) obs_kind = RADIOSONDE_V_WIND_COMPONENT
     if(obstype == 230 .or. obstype == 231) obs_kind = AIRCRAFT_V_WIND_COMPONENT
     if(obstype == 233                    ) obs_kind = ACARS_V_WIND_COMPONENT
     if(obstype == 242 .or. obstype == 243) obs_kind = SAT_V_WIND_COMPONENT
     if(obstype == 245 .or. obstype == 246) obs_kind = SAT_V_WIND_COMPONENT
     if(obstype == 252 .or. obstype == 253) obs_kind = SAT_V_WIND_COMPONENT
     if(obstype == 255                    ) obs_kind = SAT_V_WIND_COMPONENT
     if(obstype == 280 .or. obstype == 282) obs_kind = MARINE_SFC_V_WIND_COMPONENT
     if(obstype == 281 .or. obstype == 284) obs_kind = LAND_SFC_V_WIND_COMPONENT
   endif

   if (obs_kind < 0) then
      ! the "real" fix if the record type is not found might actually be to
      ! accept all record types within valid ranges, and depend on the first
      ! preprocessing steps (in the prepbufr converter) to remove obs record
      ! types which are not desired.  for now, avoid giving them the wrong type
      ! and quietly loop.
      if (obs_prof == 5) then
         if (debug) write(*,*) 'unwanted moisture obs_prof, skipping', obs_prof, obstype, zob2
         iskip(fail_moisttype) = iskip(fail_moisttype) + 1
      else
         if (debug) write(*,*) 'unrecognized obs_prof or obstype, skipping', obs_prof, obstype
         iskip(fail_badkind) = iskip(fail_badkind) + 1
      endif
      cycle obsloop 
   endif

!   check to see if this observation is desired
!------------------------------------------------------------------------------

   ! if select_obs is true, we are going to include all observations
   ! and we skip the selection code below.
   if(select_obs) then

      ! assume we are going to ignore this obs, unless it is
      ! specifically included by one of the selections below.
      pass = .true.

      ! select the specific NCEP obs types
      if( (ADPUPA .and. (subset =='ADPUPA')) .or. &
          (AIRCAR .and. (subset =='AIRCAR')) .or. &
          (AIRCFT .and. (subset =='AIRCFT')) .or. &
          (SATEMP .and. (subset =='SATEMP')) .or. &
          (SFCSHP .and. (subset =='SFCSHP')) .or. &
          (ADPSFC .and. (subset =='ADPSFC')) .or. &
          (SATWND .and. (subset =='SATWND'))       ) then

         ! then select the obs kind requested
         if( (obs_T                 .and. (obs_kind_gen == QTY_TEMPERATURE ))      .or. &
             (obs_U                 .and. (obs_kind_gen == QTY_U_WIND_COMPONENT )) .or. &
             (obs_V                 .and. (obs_kind_gen == QTY_V_WIND_COMPONENT )) .or. &
             (obs_PS                .and. (obs_kind_gen == QTY_SURFACE_PRESSURE))  .or. &
             (obs_QV                .and. (obs_kind_gen == QTY_SPECIFIC_HUMIDITY)) .or. &
             (include_relative_humidity .and. (obs_kind_gen == QTY_RELATIVE_HUMIDITY)) .or. &
             (include_dewpoint          .and. (obs_kind_gen == QTY_DEWPOINT)) ) then
             pass = .false.
         endif

      endif

      ! if pass is still true, we want to ignore this obs.
      if(pass) then
         if (debug) write(*,*) 'obs skipped because not on wanted list.  subset, obs_kind = ', subset, obs_kind_gen
         iskip(fail_notwanted) = iskip(fail_notwanted) + 1
         cycle obsloop 
      endif

   endif

!   process this observation
!------------------------------------------------------------------------------

   obs_num = obs_num + 1

   ! print a reassuring message after every Nth processed obs.
   ! if requested, print in the form of a timestamp.  
   ! the default is just a plain string with the current obs count.
   if(mod(obs_num, print_every_Nth) == 0) then
       write(label, *) 'obs count = ', obs_num
       if (print_timestamps) then
          call timestamp(string1=label, pos='brief')
       else
          write(*,*) trim(label)
       endif
   endif
   if(obs_num == max_num_obs) then
      print*, 'Max limit for observation count reached.  Increase value in namelist'
      stop
   endif

   ! set vertical coordinate for upper-air observations
   if (subset == 'AIRCAR' .or. subset == 'AIRCFT' .or. &
       subset == 'SATEMP' .or.                         &
       subset == 'SATWND' .or. subset == 'ADPUPA' ) then
       vloc = lev*100.0_r8          ! convert from mb to Pascal
       which_vert = VERTISPRESSURE
   endif

   ! set vertical coordinate for surface observations
   if (subset == 'ADPSFC' .or. subset == 'SFCSHP') then
     vloc = lev
     ! some obs have elevation of 1e12 - toss those.
     if ( subset == 'SFCSHP' .and. vloc > 4000.0_r8) then
        iskip(fail_badvert) = iskip(fail_badvert) + 1
        cycle obsloop
     endif
     which_vert = VERTISSURFACE
   endif

   ! set obs value and error if necessary
   if ( obs_kind == LAND_SFC_ALTIMETER .or. obs_kind == MARINE_SFC_ALTIMETER &
        .or. obs_kind == RADIOSONDE_SURFACE_ALTIMETER ) then
      vloc = lev                  ! station height, not used now for Ps obs
      which_vert = VERTISSURFACE
      obs_value  = compute_altimeter(zob, vloc)  !  altimeter is hPa
   elseif(obs_kind_gen == QTY_SURFACE_PRESSURE) then
      obs_value = zob * 100.0_r8  !  for Ps variable only in Pascal
      vloc = lev                  ! station height, not used now for Ps obs
      which_vert = VERTISSURFACE
   else if(obs_kind_gen == QTY_SPECIFIC_HUMIDITY) then
      obs_err = obs_err*1.0e-3_r8
      obs_value = zob*1.0e-3_r8     !  for Q variable to kg/kg
   else
      obs_value = zob               !  for T, U, V, RH, Tdew
   endif

   if (time < 24.0_r8) then
      if (wrap) then
        rday = gday - 1
      else
        rday = gday
      endif
      seconds = time * 3600
   else
      rday = gday 
      seconds = (time - 24.0_r8) * 3600
   endif

   ! fill the time_obs variable to pass into add routine
   time_obs = set_time(seconds, rday)

   !   create the obs_def for this observation, add to sequence

   call create_3d_obs(lat, lon, vloc, which_vert, obs_value, &
                      obs_kind, obs_err, rday, seconds, real(iqc,r8), obs)
   call add_obs_to_seq(seq, obs, time_obs, prev_obs, prev_time, first_obs)

end do obsloop

200 continue

print*, 'num obs used = ', obs_num, ' total obs skipped = ', sum(iskip)
do i=1, num_fail_kinds
  if (iskip(i) >  0) print *, iskip(i), 'skipped because ', skip_reasons(i)
enddo

end subroutine construct_obs_sequence


end program prepbufr_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

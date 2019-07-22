! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program wod_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     wod_to_obs - program that reads a list of World Ocean Datatbase obs
!                  profiles of ocean temperature and salinity in packed
!                  ascii format and writes a DART observation sequence file. 
!
!     created 19-oct-2009, based on the GTSPP reader.   nsc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                             increment_time, get_time, set_date, operator(-),  &
                             print_date, operator(+), leap_year, operator(>)
use    utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                             check_namelist_read, nmlfileunit, do_output,    &
                             get_next_filename, error_handler, E_ERR, E_MSG, &
                             find_textfile_dims, finalize_utilities,         &
                             open_file, close_file
use     location_mod, only : VERTISHEIGHT, set_location
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,       &
                             static_init_obs_sequence, init_obs, destroy_obs, &
                             write_obs_seq, init_obs_sequence, get_num_obs,   &
                             insert_obs_in_seq, destroy_obs_sequence,         &
                             set_copy_meta_data, set_qc_meta_data, set_qc,    & 
                             set_obs_values, set_obs_def, insert_obs_in_seq
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location, &
                             set_obs_def_key
use obs_kind_mod,     only : get_name_for_type_of_obs

! FIXME: i don't have all the actual kinds yet (bottle? underway? )
! but it's closer than before.  must have obs_def_ocean_mod.f90 in preprocess list.

use     obs_kind_mod, only : &
         SALINITY,                      QTY_SALINITY,              &
         TEMPERATURE,                   QTY_TEMPERATURE,           &
         U_CURRENT_COMPONENT,           QTY_U_CURRENT_COMPONENT,   &
         V_CURRENT_COMPONENT,           QTY_V_CURRENT_COMPONENT,   &
         SEA_SURFACE_HEIGHT,            QTY_SEA_SURFACE_HEIGHT,    &
         SEA_SURFACE_PRESSURE,          QTY_SEA_SURFACE_PRESSURE,  &
         ARGO_U_CURRENT_COMPONENT,      ARGO_V_CURRENT_COMPONENT,   &
         ARGO_SALINITY,                 ARGO_TEMPERATURE,           &
         ADCP_U_CURRENT_COMPONENT,      ADCP_V_CURRENT_COMPONENT,   &
         ADCP_SALINITY,                 ADCP_TEMPERATURE,           &
         FLOAT_SALINITY,                FLOAT_TEMPERATURE

use obs_kind_mod, only : &
         DRIFTER_U_CURRENT_COMPONENT,   DRIFTER_V_CURRENT_COMPONENT,   &
         DRIFTER_SALINITY,              DRIFTER_TEMPERATURE,           &
         GLIDER_U_CURRENT_COMPONENT,    GLIDER_V_CURRENT_COMPONENT,    &
         GLIDER_SALINITY,               GLIDER_TEMPERATURE,            &
         MOORING_U_CURRENT_COMPONENT,   MOORING_V_CURRENT_COMPONENT,   &
         MOORING_SALINITY,              MOORING_TEMPERATURE,           &
         MOORING_PRESSURE

use obs_kind_mod, only : &
         CTD_SALINITY,                  CTD_TEMPERATURE,               &
         TCTD_SALINITY,                 TCTD_TEMPERATURE,              &
         XCTD_SALINITY,                 XCTD_TEMPERATURE,              &
         STD_SALINITY,                  STD_TEMPERATURE,               &
         MBT_TEMPERATURE,               XBT_TEMPERATURE,               &
         DBT_TEMPERATURE,               APB_TEMPERATURE,               &
         MBT_SALINITY,                  XBT_SALINITY,                  &
         DBT_SALINITY,                  APB_SALINITY,                  &
         BOTTLE_TEMPERATURE,            BOTTLE_SALINITY,               &
         DOPPLER_U_CURRENT_COMPONENT,   DOPPLER_V_CURRENT_COMPONENT,   &
         DOPPLER_W_CURRENT_COMPONENT,   QTY_W_CURRENT_COMPONENT,      &
         SATELLITE_MICROWAVE_SST,       SATELLITE_INFRARED_SST,        &
         SATELLITE_SSH,                 SATELLITE_SSS,                 &
         HFRADAR_RADIAL_VELOCITY,       QTY_VELOCITY


! not clean interface; going for working code first.
use WOD_read_routines_mod, only : WODREADDART, depth, temp, ierror, iderror,  &
                      maxlevel, maxcalc, kdim, maxtcode, maxtax, maxpinf,     &
                      bmiss, isec, sechead



implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=129) :: msgstring, next_infile, cdummy
integer :: j, k, nfiles, num_new_obs, castid, l
integer :: oday, osec                   
integer :: obsyear, obsmonth, obsday, obssec
integer :: obs_num, io, iunit, filenum, dummy
integer :: funit, levels, istdlev, nvar, nsecond, ieof
integer :: ip2(0:maxlevel), cast, itype
logical :: file_exist, did_obs, from_list = .false.
logical :: have_temp, have_salt
real(r8) :: terr, serr,         & 
            phs, obs_val(1),  &
            d_qc(1)
real :: dtime, lato, lono
real(r8) :: obslat, obslon, obsdepth

! platform codes - reshape required by latest XLF compiler.
integer :: ptype(2, 16) = reshape( (/  &
   MBT_TEMPERATURE,     MBT_SALINITY,      &   ! ptype(1)  = MBT
   XBT_TEMPERATURE,     XBT_SALINITY,      &   ! ptype(2)  = XBT
   DBT_TEMPERATURE,     DBT_SALINITY,      &   ! ptype(3)  = DBT
   CTD_TEMPERATURE,     CTD_SALINITY,      &   ! ptype(4)  = CTD
   STD_TEMPERATURE,     STD_SALINITY,      &   ! ptype(5)  = STD
   XCTD_TEMPERATURE,    XCTD_SALINITY,     &   ! ptype(6)  = XCTD
   BOTTLE_TEMPERATURE,  BOTTLE_SALINITY,   &   ! ptype(7)  = bottle
   0,                   0,                 &   ! ptype(8)  = underway
   FLOAT_TEMPERATURE,   FLOAT_SALINITY,    &   ! ptype(9)  = profiling float
   MOORING_TEMPERATURE, MOORING_SALINITY,  &   ! ptype(10) = moored buoy
   DRIFTER_TEMPERATURE, DRIFTER_SALINITY,  &   ! ptype(11) = drifting buoy
   TCTD_TEMPERATURE,    TCTD_SALINITY,     &   ! ptype(12) = towed CTD
   APB_TEMPERATURE,     APB_SALINITY,      &   ! ptype(13) = animal
   0,                   0,                 &   ! ptype(14) = bucket
   GLIDER_TEMPERATURE,  GLIDER_SALINITY,   &   ! ptype(15) = glider
   MBT_TEMPERATURE,     MBT_SALINITY  /),  &   ! ptype(16) = microBT
   (/ 2, 16 /)  )  ! reshape 1d array into (2,16)


type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, delta_time, prev_time

integer :: histbin(-1:26) = 0  ! for time debug
integer :: histcount = 0       ! for time debug
integer :: temp_type, salt_type, good_temp, good_salt, bad_temp, bad_salt
integer :: salt_qc(10), temp_qc(10)

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=128) :: wod_input_file     = 'XBTS2005'
character(len=128) :: wod_input_filelist = ''
character(len=128) :: wod_out_file       = 'obs_seq.wod'
integer            :: avg_obs_per_file   = 500000
logical            :: debug              = .false.
logical            :: timedebug          = .false.
logical            :: print_qc_summary   = .true.
integer            :: max_casts          = -1
logical            :: no_output_file     = .false.
integer            :: print_every_nth_cast = -1
real(r8)           :: temperature_error  = 0.5  ! degrees C
real(r8)           :: salinity_error     = 0.5  ! g/kg
integer            :: start_month        = 1
integer            :: end_month          = 12

namelist /wod_to_obs_nml/ &
   wod_input_file, wod_input_filelist, wod_out_file,   &
   avg_obs_per_file, debug, max_casts, no_output_file, &
   print_every_nth_cast, print_qc_summary,             &
   temperature_error, salinity_error, timedebug,       &
   start_month, end_month

! start of executable code

obs_num = 1
d_qc(1) = 0.0_r8

! time in files is year, month, day, fraction of day.
! day is sometimes = 0 (? not sure, am setting it to day 1).
! time often missing; am setting that to 0Z with less
! trepidation.
call set_calendar_type(GREGORIAN)
prev_time = set_date(4000, 1, 1)   ! must be something later than all obs

!  read the necessary parameters from input.nml
call initialize_utilities()
call find_namelist_in_file("input.nml", "wod_to_obs_nml", iunit)
read(iunit, nml = wod_to_obs_nml, iostat = io)

! Record the namelist values used for the run 
if (do_output()) write(nmlfileunit, nml=wod_to_obs_nml)

! any needed namelist checks for sanity:

! cannot have both a single filename and a list; the namelist must
! shut one off.
if (wod_input_file /= '' .and. wod_input_filelist /= '') then
  call error_handler(E_ERR, 'wod_to_obs',                     &
                     'One of wod_input_file or filelist must be NULL', &
                     source, revision, revdate)
endif
if (wod_input_filelist /= '') from_list = .true.

! need to know a reasonable max number of obs that could be added here.
if (from_list) then
   call find_textfile_dims(wod_input_filelist, nfiles, dummy)
   num_new_obs = avg_obs_per_file * nfiles
else
   num_new_obs = avg_obs_per_file
endif

! FIXME: is this a good idea?  to append to an existing file?
!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
inquire(file=wod_out_file, exist=file_exist)
if (file_exist) then
   print *, "replacing existing obs_seq file ", trim(wod_out_file)
else
   print *, "creating obs_seq file ", trim(wod_out_file)
endif
print *, "max entries = ", num_new_obs
call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
do k = 1, num_copies
   call set_copy_meta_data(obs_seq, k, 'WOD observation')
end do
do k = 1, num_qc
   call set_qc_meta_data(obs_seq, k, 'WOD QC')
end do

did_obs = .false.

! main loop that does either a single file or a list of files

filenum = 1
fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(wod_input_filelist, filenum)
   else
      next_infile = wod_input_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop
  
   have_temp = .false.
   have_salt = .false.
   !FIXME: more types?  u, v, w  current, ?

   !  open the next profile file
   funit = open_file(next_infile, action='READ')

print *, ' '
print *, 'opening file ', trim(next_infile)

   ! reset anything that's per-file
   good_temp = 0
   bad_temp = 0
   good_salt = 0
   bad_salt = 0
   temp_qc(:) = 0
   salt_qc(:) = 0

   cast = 1
   castloop: do    ! until out of data
   call WODreadDART(funit,obsyear,obsmonth,obsday, &
              dtime,lato,lono,levels,istdlev,nvar,ip2,nsecond, &
              bmiss,castid,ieof)

!if (ieof /= 0) print *, 'ieof is ', ieof

   ! this comes back as 1 when the file has all been read.
   if (ieof == 1) exit castloop

   ! see what data we have
   if (nvar == 0) then
      goto 20 ! inc counter, cycle castloop
   endif

!print *, 'obsyear, obsmonth, obsday = ', obsyear, obsmonth, obsday
!print *, 'levels = ', levels
!print *, 'lato, lono = ', lato, lono

   ! need to gather:
   !   time (year, month, day, time), lat, lon, depth, obs type, ?.

   ! start out with converting the real time.
   ! do error check first.  if fail, cycle here
   ! this routine alters bad days and times to be valid, or
   ! returns a failure and we loop here.
   if (.not. date_ok(obsyear, obsmonth, obsday, &
                     dtime, bmiss, castid)) goto 20 
 
   ! convert fractional hours (0-24) to integer days and seconds
   !  and add on to reference time.  dtime was set to 0 if it was
   !  missing in the date_ok() routine.
   obssec = int(dtime * 3600.0)
   delta_time = set_time(obssec)
   obs_time = set_date(obsyear, obsmonth, obsday) + delta_time
   call get_time(obs_time,  osec, oday)
!print *, 'oday, osec = ', oday, osec
!call print_date(obs_time)
   
   if (debug) then
      print *, ' --- '
      print *, 'cast number: ', castid
      print *, 'obsyear,month,day,time = ', obsyear, obsmonth, obsday, dtime
      print *, 'lato, lono = ', lato, lono
   endif

   have_temp = .false.
   have_salt = .false.
   do k=1, nvar
      if (ip2(k) == 1) have_temp = .true.
      if (ip2(k) == 2) have_salt = .true.
   enddo

   temp_type = 0
   salt_type = 0
   do k=1, nsecond
      if (isec(k) == 29) then
         itype = int(sechead(29))
         if (itype > 0 .and. itype <= 16) then
            temp_type = ptype(1, itype)
            salt_type = ptype(2, itype)
         endif
         if (itype > 0 .and. temp_type == 0) print *, 'temp=0, itype = ', itype
         if (itype > 0 .and. salt_type == 0) print *, 'salt=0, itype = ', itype
      endif
   enddo

! debug
   if (debug) then
      print *, 'num levels: ', levels
      print *, 'has temp, salinity: ', have_temp, have_salt
      if (have_temp) then
         cdummy = get_name_for_type_of_obs(temp_type)
         print *, 'temp type/name: ', temp_type, trim(cdummy)
      endif 
      if (have_salt) then
         cdummy = get_name_for_type_of_obs(salt_type)
         print *, 'salt type/name: ', salt_type, trim(cdummy)
      endif 
      call print_date(obs_time, 'obs time')
   endif
   
 ! FIXME: these are best-guess only.
   terr = temperature_error            ! degrees C
   serr = salinity_error / 1000.0_r8   ! g/kg, convert to kg/kg

   if (lono < 0.0_r8) lono = lono + 360.0_r8
   obslon = lono 
   obslat = lato
   if((obslon < 0.0_r8 .or. obslon > 360.0_r8) .or. (obslat < -90.0_r8 .or. obslat > 90.0_r8)) then
      print *, 'FSC: longitude (',obslon,') or latitude (',obslat,') is not within range. cycle castloop'
      goto 20 ! inc counter, cycle castloop
   endif

   if (have_temp) then
      if (ierror(1) == 0) then
         good_temp = good_temp + 1
      else
         bad_temp = bad_temp + 1
         temp_qc(ierror(1)) = temp_qc(ierror(1)) + 1
      endif
   endif
   
   if (have_salt) then
      if (ierror(2) == 0) then
         good_salt = good_salt + 1
      else
         bad_salt = bad_salt + 1
         salt_qc(ierror(2)) = salt_qc(ierror(2)) + 1
      endif
   endif
   
   if (debug) then
      if (ierror(1) /= 0) print *, 'whole temp cast discarded, ierror == ', ierror(1)
      if (ierror(2) /= 0) print *, 'whole salt cast discarded, ierror == ', ierror(2)
   endif

   ! the incoming files are yearly collections.  allow this converter to create
   ! output files which are only partial years by only processing selected months.
   if (obsmonth < start_month .or. obsmonth > end_month) then
      goto 20 ! inc counter, cycle castloop
   endif

   obsloop: do k = 1, levels
   
     obsdepth = depth(k)

      if (debug) then
         write(*, "(A,F6.0,F8.3,F4.0,F8.3,F4.0)") "depth, temp/salt value,qc: ",&
                    depth(k), temp(k,1), iderror(k,1), temp(k, 2), iderror(k,2)
      endif

     ! ierror is whole cast error
     if (have_temp .and. ierror(1) == 0 .and. &
         temp_type > 0 .and. temp(k, 1) /= bmiss &
         .and. (.not. no_output_file)) then   

         ! set location - incoming obs are -180 to 180 in longitude;
         ! dart wants 0 to 360.  (also, r8)
         call fill_obs(obs, obs_num, obslon, obslat, obsdepth, temp_type, &
                 obs_time, terr, real(temp(k, 1), r8), real(iderror(k, 1), r8))

         call add_obs(obs_seq, obs, obs_time, prev_obs, prev_time)

         if (.not. did_obs) did_obs = .true.
      endif


      ! ierror is whole cast
      if (have_salt .and. ierror(2) == 0 .and. &
          salt_type > 0 .and. temp(k, 2) /= bmiss &
         .and. (.not. no_output_file)) then   

         call fill_obs(obs, obs_num, obslon, obslat, obsdepth, salt_type, &
               obs_time, serr, temp(k, 2) / 1000.0_r8, real(iderror(k, 2), r8))

         call add_obs(obs_seq, obs, obs_time, prev_obs, prev_time)

         if (.not. did_obs) did_obs = .true.
      endif 

   end do obsloop

 20 continue
  if (print_every_nth_cast > 0) then
      if (mod(cast, print_every_nth_cast) == 0) then
         print *, 'processing cast ', cast
      endif
  endif

  cast = cast + 1
  if (max_casts > 0 .and. cast >= max_casts) exit castloop

  end do castloop

  call close_file(funit)
  filenum = filenum + 1

  if (print_qc_summary) then
     call error_handler(E_MSG, '', '')
     write(msgstring, *) 'input data filename: ', trim(next_infile)
     call error_handler(E_MSG, '', msgstring)
     call error_handler(E_MSG, '', '')
     write(msgstring, *) 'total casts, total profiles: ', cast, &
                          good_temp+bad_temp+good_salt+bad_salt
     call error_handler(E_MSG, '', msgstring)
     write(msgstring, *) 'total/good/bad temperature profiles: ', &
                          good_temp+bad_temp, good_temp, bad_temp
     call error_handler(E_MSG, '', msgstring)
     write(msgstring, *) 'total/good/bad salinity    profiles: ', &
                          good_salt+bad_salt, good_salt, bad_salt
     call error_handler(E_MSG, '', msgstring)
     call error_handler(E_MSG, '', '')
      do j=1, 10
         if (temp_qc(j) > 0) then
            write(msgstring, *) 'temp qc of ', j, ' found ', temp_qc(j), ' times'
            call error_handler(E_MSG, '', msgstring)
         endif
      enddo
      call error_handler(E_MSG, '', '')
      do j=1, 10
         if (salt_qc(j) > 0) then
            write(msgstring, *) 'salt qc of ', j, ' found ', salt_qc(j), ' times'
            call error_handler(E_MSG, '', msgstring)
         endif
      enddo
      call error_handler(E_MSG, '', '')
   endif

end do fileloop

! done with main loop.  if we added any obs to the sequence, write it out.
if (did_obs) then
print *, 'ready to write, nobs = ', get_num_obs(obs_seq)
   if (get_num_obs(obs_seq) > 0) &
      call write_obs_seq(obs_seq, wod_out_file)

   ! minor stab at cleanup, in the off chance this will someday get turned
   ! into a subroutine in a module.  probably not all that needs to be done,
   ! but a start.
   call destroy_obs(obs)
   !call destroy_obs(prev_obs)   ! is this identical to obs?
   ! get core dumps here, not sure why?
   if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
endif

! histogram of times, debug only
if (timedebug) then
   print *, "total obs examined = ", histcount
   print *, "bin -1 is < 0"
   print *, "bin 0 is == 0"
   print *, "bin 25 is > 24 but != 99.9"
   print *, "bin 26 is == 99.9"
   print *, " "
   print *, 'bin    count'
   do l=lbound(histbin, 1), ubound(histbin, 1)
      print *, l, histbin(l)
   enddo
endif

call error_handler(E_MSG,'wod_to_obs','Finished successfully.',source,revision,revdate)
call finalize_utilities()


! END OF MAIN ROUTINE

contains

! local subroutines/functions follow

! subroutine to fill an obs.  assumes loc is 3d, and vert is in meters.
subroutine fill_obs(obs, onum, olon, olat, odepth, otype, otime, oerr, &
                    oval, oqc)
 type(obs_type),  intent(inout) :: obs
 integer,         intent(inout) :: onum
 real(r8),        intent(in)    :: olon, olat, odepth, oerr, oval, oqc
 integer,         intent(in)    :: otype
 type(time_type), intent(in)    :: otime

type(obs_def_type) :: obs_def
real(r8)           :: valarr(1), qcarr(1)

call set_obs_def_location(obs_def, &
                          set_location(olon, olat, odepth, VERTISHEIGHT))
call set_obs_def_type_of_obs(obs_def, otype)
call set_obs_def_time(obs_def, otime)
    
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def_key(obs_def, onum)
call set_obs_def(obs, obs_def)

valarr(1) = oval
call set_obs_values(obs, valarr)
qcarr(1)  = oqc
call set_qc(obs, qcarr)

onum = onum + 1

end subroutine fill_obs
    
! add an obs to the sequence.  if prev_time same or earlier
! than obs time, insert with search starting from prev obs.
subroutine add_obs(seq, obs, obs_time, prev_obs, prev_time)

! FIXME ... 'seq' in argument list never used ... 
 type(obs_sequence_type),   intent(inout) :: seq
 type(obs_type),            intent(inout) :: obs, prev_obs
 type(time_type),           intent(in)    :: obs_time
 type(time_type),           intent(inout) :: prev_time

logical, save :: first_obs = .true.

if (first_obs .or. prev_time > obs_time) then
   call insert_obs_in_seq(obs_seq, obs)
   first_obs = .false.
else
   call insert_obs_in_seq(obs_seq, obs, prev_obs)
endif

prev_obs  = obs
prev_time = obs_time

end subroutine add_obs

! date check, since a lot of these obs seem to have bad times/dates.
! modify testmonth, testday, dtime if we know what the right answer is.
! return .TRUE. if date ok, .FALSE. if bad and not fixable.
function date_ok(testyear, testmonth, testday, dtime, bmiss, castid)
 integer,  intent(inout) :: testyear, testmonth, testday
 real,     intent(inout) :: dtime    ! note real, not real(r8), to match 
 real,     intent(in)    :: bmiss    ! WOD-supplied code which uses reals
 integer,  intent(in)    :: castid
 
 logical :: date_ok

integer :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
type(time_type) :: test_time

! at least 2 files have dates of 2/29 on non-leap years.  test for
! that and fix those dates.  also found a 9/31, apparently without
! a bad cast qc.   the 2/29s were examined and determined to be a bad
! conversion - should be mar 1  the 9/31s were recording errors and
! based on the other obs, should have been 9/30.

10 format (A,I9,A,I6,I4,I4,F6.2)

! bad month - return fail
if (testmonth <= 0 .or. testmonth > 12) then
   write(*, 10) 'cast number: ', castid, ' date, dtime: ', testyear, testmonth, testday, dtime
   write(*, *)  ' bad month number, discarding observation'
   date_ok = .FALSE.
   return 
endif

! try to estimate the time shifting error impacts from old code.
if (timedebug) then
   histcount = histcount + 1
   if (dtime == bmiss) then
      ! no count.  it is ok since we discarded it before
   else if (dtime < 0.0) then
      histbin(-1) = histbin(-1) + 1
   else if (dtime == 0.0) then
      histbin(0) = histbin(0) + 1
   else if (dtime == 99.9) then
      histbin(26) = histbin(26) + 1
   else if (dtime > 24.0) then
      histbin(25) = histbin(25) + 1
   else
      histbin(int(dtime)) = histbin(int(dtime)) + 1
   endif
endif

! check time - if standard missing value, don't squawk.
! otherwise, note the cast number and original value, and
! then set to 0, or fail.
if (dtime == bmiss) then
   dtime = 0.0
else if (dtime == 99.90) then
   ! seems to be an alternative missing value in some few cases.
   ! got confirmation from WOD providers that this was an incorrect
   ! missing value and future versions of the WOD will use the standard
   ! 'missing_value' read out of the file (bmiss in this subroutine, 
   ! appears to be something like -999.99)
   write(*, 10) 'cast number: ', castid, ' date, dtime: ', testyear, testmonth, testday, dtime
   write(*, *)  ' setting to 0 and continuing'
   dtime = 0.0
else if (dtime < 0) then
   write(*, 10) 'cast number: ', castid, ' date, dtime: ', testyear, testmonth, testday, dtime
   write(*, *)  ' setting to 0 and continuing'
   dtime = 0.0
else if (dtime > 24.0) then
   write(*, 10) 'cast number: ', castid, ' date, dtime: ', testyear, testmonth, testday, dtime
   ! some obs seem to have more than 24 hours in the time variable.
   ! got confirmation from WOD providers that some appear to be keypunch
   ! entry errors.  this field should never be > 24.  for now, set to 0
   ! because the date is believed to be good.
   !write(*, *)  ' discarding obs'
   !date_ok = .FALSE.
   !return 
   write(*, *)  ' setting to 0 and continuing'
   dtime = 0.0
endif

! good day - return success
if (testday > 0 .and. testday <= days_per_month(testmonth)) then
   date_ok = .TRUE.
   return 
endif
 
! don't do the year test unless this is leap day.
if ((testmonth == 2) .and. (testday == 29)) then
   ! jan 1st always exists.  set this so we can test for leap year.
   test_time = set_date(testyear, 1, 1)
   if (.not. leap_year(test_time)) then
      write(*, 10) 'cast number: ', castid, ' date, dtime: ', testyear, testmonth, testday, dtime
      write(*, *)  ' does not exist in this non-leap year; setting to mar 1'
      ! convert this to mar 1 (date conversion error)
      testmonth = 3
      testday = 1
   endif
   date_ok = .TRUE.
   return 
endif

! day 0 is possibly a processing error, but actual day is currently unknown.
! i have been advised to discard these, which i am reluctantly going to do.
if (testday == 0) then
   write(*, 10) 'cast number: ', castid, ' date, dtime: ', testyear, testmonth, testday, dtime
   write(*, *)  ' does not exist; discarding obs'
   date_ok = .FALSE.
   return 
   ! alternative code - keep obs and set day to 1 and time to 0.
   !write(*, *)  ' setting day to 1, dtime to 0'
   !testday = 1
   !dtime = 0.0
   !date_ok = .TRUE.
   return 
endif

! if only 1 > end of month, set to end of month and succeed.
! if more than 1 > end of month, fail.
! found other files with 4/31 and 6/39 as dates. 
if (testday > days_per_month(testmonth)+1) then
   write(*, 10) 'cast number: ', castid, ' date, dtime: ', testyear, testmonth, testday, dtime
   write(*, *)  ' does not exist; discarding obs'
   date_ok = .FALSE.
   return 
endif

! not sure this is a good idea; maybe should discard these obs?
if (testday == days_per_month(testmonth)+1) then
   write(*, 10) 'cast number: ', castid, ' date, dtime: ', testyear, testmonth, testday, dtime
   write(*, *)  ' does not exist; assuming rollover to day 1, month+1, dtime to 0'
   dtime = 0.0
   testday = 1
   testmonth = testmonth + 1
   if (testmonth > 12) then
      testmonth = 1
      testyear = testyear + 1
   endif
   date_ok = .TRUE.
   return 
endif

! shouldn't get here.  assume bad date.
write(*, 10) 'cast number: ', castid, ' date, dtime: ', testyear, testmonth, testday, dtime
write(*, *)  ' bad date; discarding obs'
date_ok = .FALSE.
return 

end function date_ok

end program wod_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

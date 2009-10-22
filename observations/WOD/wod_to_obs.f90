! Data Assimilation Research Testbed -- DART
! Copyright 2004-2008, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     wod_to_obs - program that reads a list of World Ocean Datatbase obs
!                  profiles of ocean temperature and salinity in packed
!                  ascii format and writes a DART observation sequence file. 
!
!     created 19-oct-2009, based on the GTSPP reader.   nsc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program wod_to_obs

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                             increment_time, get_time, set_date, operator(-),  &
                             print_date, operator(+)
use    utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                             check_namelist_read, nmlfileunit, do_output,    &
                             get_next_filename, error_handler, E_ERR, E_MSG, &
                             nc_check, find_textfile_dims, finalize_utilities, &
                             timestamp, open_file, close_file
use     location_mod, only : VERTISHEIGHT, set_location
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,       &
                             static_init_obs_sequence, init_obs, destroy_obs, &
                             write_obs_seq, init_obs_sequence, get_num_obs,   &
                             insert_obs_in_seq, destroy_obs_sequence,         &
                             set_copy_meta_data, set_qc_meta_data, set_qc,    & 
                             set_obs_values, set_obs_def, insert_obs_in_seq
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location, &
                             set_obs_def_key

! FIXME: i don't have all the actual kinds yet (bottle? underway? )
! but it's closer than before.  must have obs_def_ocean_mod.f90 in preprocess list.

use     obs_kind_mod, only : &
         SALINITY,                      KIND_SALINITY,              &
         TEMPERATURE,                   KIND_TEMPERATURE,           &
         U_CURRENT_COMPONENT,           KIND_U_CURRENT_COMPONENT,   &
         V_CURRENT_COMPONENT,           KIND_V_CURRENT_COMPONENT,   &
         SEA_SURFACE_HEIGHT,            KIND_SEA_SURFACE_HEIGHT,    &
         SEA_SURFACE_PRESSURE,          KIND_SEA_SURFACE_PRESSURE,  &
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
         DOPPLER_W_CURRENT_COMPONENT,   KIND_W_CURRENT_COMPONENT,      &
         SATELLITE_MICROWAVE_SST,       SATELLITE_INFRARED_SST,        &
         SATELLITE_SSH,                 SATELLITE_SSS,                 &
         CODAR_RADIAL_VELOCITY,         KIND_VELOCITY


! not clean interface; going for working code first.
use WOD_read_routines_mod, only : WODREADDART, depth, temp, ierror, iderror,  &
                      maxlevel, maxcalc, kdim, maxtcode, maxtax, maxpinf, bmiss, &
                      isec, sechead



implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"


integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=129) :: meta_data, msgstring, next_infile
character (len=80)  :: name
character (len=19)  :: datestr
character (len=6)   :: subset
integer :: rcode, ncid, varid, k, nfiles, num_new_obs
integer :: aday, asec, dday, dsec, oday, osec                   
integer :: obsyear, obsmonth, obsday, obshour, obsmin, obssec
integer :: zloc, obs_num, io, iunit, filenum, dummy, i_qc, nc_rc
integer :: funit, levels, istdlev, nvar, nsecond, ieof
integer :: ip2(0:maxlevel), cast, itype
logical :: file_exist, first_obs, did_obs, from_list = .false.
logical :: have_temp, have_salt
real(r8) :: hght_miss, refr_miss, azim_miss, terr, serr,         & 
            qc, hghto, refro, azimo, wght, nx, ny,   & 
            nz, ds, htop, rfict, obsval, phs, obs_val(1), qc_val(1),  &
            glat, glon, d_qc(1)
real :: dtime, lato, lono
real(r8) :: obslat, obslon, obsdepth

! platform codes
integer :: ptype(2, 16) = (/  &
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
   MBT_TEMPERATURE,     MBT_SALINITY  /)       ! ptype(16) = microBT


real(r8), allocatable :: lat(:), lon(:), dep(:), err(:) !, d_qc(:)

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, delta_time

integer :: temp_type, salt_type, bad_casts

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=128) :: wod_input_file     = 'XBTS2005'
character(len=128) :: wod_input_filelist = ''
character(len=128) :: wod_out_file       = 'obs_seq.wod'
integer            :: avg_obs_per_file   = 500000
logical            :: debug              = .false.
integer            :: max_casts          = -1
integer            :: print_every_nth_cast = -1

namelist /wod_to_obs_nml/ wod_input_file,                   &
                          wod_input_filelist, wod_out_file, &
                          avg_obs_per_file, debug, max_casts, &
                          print_every_nth_cast

! start of executable code

obs_num = 1
d_qc(1) = 0.0_r8

! time in files is year, month, day, fraction of day.
! day is sometimes = 0 (? not sure, am setting it to day 1).
! time often missing; am setting that to 0Z with less
! trepidation.
call set_calendar_type(GREGORIAN)

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
   meta_data = 'WOD observation'
   call set_copy_meta_data(obs_seq, k, meta_data)
end do
do k = 1, num_qc
   meta_data = 'WOD QC'
   call set_qc_meta_data(obs_seq, k, meta_data)
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
   bad_casts = 0

   cast = 1
   castloop: do    ! until out of data
   call WODreadDART(funit,obsyear,obsmonth,obsday, &
              dtime,lato,lono,levels,istdlev,nvar,ip2,nsecond, &
              bmiss,ieof)

!if (ieof /= 0) print *, 'ieof is ', ieof

   ! this comes back as 1 when the file has all been read.
   if (ieof == 1) exit castloop

   ! if the whole cast is marked with a bad qc, loop here as well.
   ! i'm seeing a few bad dates, e.g. 2/29 on non-leap years,
   ! sept 31.  try testing the cast qc before decoding the date,
   ! since it is a fatal error in our set_time() routine to set
   ! an invalid date.
   if (ierror(1) /= 0) then
      bad_casts = bad_casts + 1
      goto 20 ! inc counter, cycle castloop
   endif

!print *, 'obsyear, obsmonth, obsday = ', obsyear, obsmonth, obsday
!print *, 'levels = ', levels
!print *, 'lato, lono = ', lato, lono

   ! FIXME: need to gather:
   !   time (year, month, day, time), lat, lon, depth, obs type, ?.

   ! start out with converting the real time.
 
   ! at least 2 files have dates of 2/29 on non-leap years.  test for
   ! that and reject those obs.  also found a 9/31, apparently without
   ! a bad cast qc.
   ! FIXME:  this is the code i wanted to use, but leap_year() isn't 
   ! implemented in the current time_manager.  hack it for now.
   !if ((.not. leap_year(obsyear)) .and. (obsmonth == 2) .and. (obsday == 29)) then
   if ((obsyear /= 2000) .and. (obsmonth == 2) .and. (obsday == 29)) then
      print *, 'date says:  ', obsyear, obsmonth, obsday
      print *, 'which does not exist in this year; skipping obs'
      goto 20 ! inc counter, cycle castloop
   endif
   if ((obsmonth == 9) .and. (obsday == 31)) then
      print *, 'date says:  ', obsyear, obsmonth, obsday
      print *, 'which does not exist in any year; skipping obs'
      goto 20 ! inc counter, cycle castloop
   endif

   ! convert to integer days and seconds, and add on to reference time.
   if (dtime == bmiss) then
      obssec = 0
   else
      obssec = int(dtime * 86400.0)
   endif
   ! fixme?
   if (obsday == 0) then
      obsday = 1
   endif
   delta_time = set_time(obssec)
   obs_time = set_date(obsyear, obsmonth, obsday) + delta_time
   call get_time(obs_time,  osec, oday)
!print *, 'oday, osec = ', oday, osec
!call print_date(obs_time)
   

   ! see what data we have
   if (nvar == 0) then
      goto 20 ! inc counter, cycle castloop
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
      print *, 'temp, salinity type: ', temp_type, salt_type
      call print_date(obs_time, 'obs time')
   endif
   
 ! FIXME: these are best-guess only.
   terr = 1.0_r8                ! temp error = 1 degrees C
   serr = 0.5_r8 / 1000.0_r8    ! salinity error = 0.5 g/kg, convert to kg/kg

   if (lono < 0.0_r8) lono = lono + 360.0_r8
   obslon = lono 
   obslat = lato

   first_obs = .true.
   
   obsloop: do k = 1, levels
   
     obsdepth = depth(k)

     ! ierror is whole cast error
     if (have_temp .and. ierror(1) == 0 .and. &
         temp_type > 0 .and. temp(k, 1) /= bmiss) then   

         ! set location - incoming obs are -180 to 180 in longitude;
         ! dart wants 0 to 360.  (also, r8)
         call set_obs_def_location(obs_def, &
                           set_location(obslon, obslat, obsdepth,VERTISHEIGHT))
         call set_obs_def_kind(obs_def, temp_type)
         call set_obs_def_time(obs_def, set_time(osec, oday))
    
         call set_obs_def_error_variance(obs_def, terr * terr)
         call set_obs_def_key(obs_def, obs_num)
         call set_obs_def(obs, obs_def)
   
!print *, 'temp = ', temp(k, 1)
         obs_val(1) = temp(k, 1)
         call set_obs_values(obs, obs_val)
         qc_val(1)  = iderror(k, 1)  ! iderror is per-depth error
         call set_qc(obs, qc_val)
    
         ! first one, insert with no prev.  otherwise, since all times will be the
         ! same for this column, insert with the prev obs as the starting point.
         ! (the first insert with no prev means it will search for the right
         ! time ordered starting point.)
         if (first_obs) then
            call insert_obs_in_seq(obs_seq, obs)
            first_obs = .false.
         else
           call insert_obs_in_seq(obs_seq, obs, prev_obs)
         endif
         obs_num = obs_num+1
         prev_obs = obs
 
         if (.not. did_obs) did_obs = .true.
      endif


      ! ierror is whole cast
      if (have_salt .and. ierror(1) == 0 .and. &
          salt_type > 0 .and. temp(k, 2) /= bmiss) then  

         call set_obs_def_location(obs_def, &
                           set_location(obslon, obslat, obsdepth,VERTISHEIGHT))
         call set_obs_def_kind(obs_def, salt_type)
         call set_obs_def_time(obs_def, set_time(osec, oday))
    
         call set_obs_def_error_variance(obs_def, serr * serr)
         call set_obs_def_key(obs_def, obs_num)
         call set_obs_def(obs, obs_def)
   
!print *, 'salt = ', temp(k, 2) / 1000.0_r8
         obs_val(1) = temp(k, 2) / 1000.0_r8
         call set_obs_values(obs, obs_val)
         qc_val(1)  = iderror(k, 2)
         call set_qc(obs, qc_val)
    
         ! first one, insert with no prev.  otherwise, since all times will be the
         ! same for this column, insert with the prev obs as the starting point.
         ! (the first insert with no prev means it will search for the right
         ! time ordered starting point.)
         if (first_obs) then
            call insert_obs_in_seq(obs_seq, obs)
            first_obs = .false.
         else
           call insert_obs_in_seq(obs_seq, obs, prev_obs)
         endif
         obs_num = obs_num+1
         prev_obs = obs
 
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

   print *, 'filename, total casts, bad casts: ', trim(next_infile), cast, bad_casts

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

call timestamp(source,revision,revdate,'end')
call finalize_utilities()

! END OF MAIN ROUTINE

! contains

! local subroutines/functions follow

! something to set the err var - make it a subroutine so we can muck with it

! subroutine to fill an obs?

end program

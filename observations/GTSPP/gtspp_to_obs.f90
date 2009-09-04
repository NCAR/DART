! Data Assimilation Research Testbed -- DART
! Copyright 2004-2008, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     gtspp_to_obs - program that reads a list of NODS GTSPP observation 
!                    profiles of ocean temperature and salinity in netcdf
!                    format and writes a DART observation sequence file. 
!
!     created 02-sept-2009, based on the GPS reader.   nsc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program gtspp_to_obs

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                             increment_time, get_time, set_date, operator(-),  &
                             print_date, operator(+)
use    utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                             check_namelist_read, nmlfileunit, do_output,    &
                             get_next_filename, error_handler, E_ERR, E_MSG, &
                             nc_check, find_textfile_dims
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
! FIXME: what actual instrument took these readings? FLOAT_xx is a placeholder
! for now.  must have obs_def_ocean_mod.f90 in the preprocess input list.
use     obs_kind_mod, only : KIND_TEMPERATURE, KIND_SALINITY, &
                             FLOAT_TEMPERATURE, FLOAT_SALINITY

use           netcdf

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
integer :: rcode, ncid, varid, ndepths, k, nfiles, num_new_obs,  &
           aday, asec, dday, dsec, oday, osec,                   &
           iyear, imonth, iday, ihour, imin, isec,               &
           zloc, obs_num, io, iunit, filenum, dummy, i_qc, nc_rc
logical :: file_exist, first_obs, did_obs, from_list = .false.
logical :: have_temp = .false., have_salt = .false.
real(r8) :: hght_miss, refr_miss, azim_miss, terr, serr,         & 
            qc, lato, lono, hghto, refro, azimo, wght, nx, ny,   & 
            nz, ds, htop, rfict, obsval, phs, obs_val(1), qc_val(1),  &
            dtime, glat, glon, d_qc(1)

real(r8), allocatable :: lat(:), lon(:), dep(:), err(:) !, d_qc(:)

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, base_time, delta_time

! initialize some values
integer, parameter :: nmaxdepths = 5000   !  max number of observation depths
real(r8) :: obs_depth(nmaxdepths)   = -1.0_r8
real(r8) :: temperature(nmaxdepths) = -888888.0_r8
real(r8) :: salinity(nmaxdepths)    = -888888.0_r8
character(len=nmaxdepths) :: t_qc = '', s_qc = ''

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=128) :: gtspp_netcdf_file     = '1234567.nc'
character(len=128) :: gtspp_netcdf_filelist = 'gtspp_to_obs_filelist'
character(len=128) :: gtspp_out_file        = 'obs_seq.gtspp'
integer            :: avg_obs_per_file      = 500

namelist /gtspp_to_obs_nml/ gtspp_netcdf_file,                     &
                            gtspp_netcdf_filelist, gtspp_out_file, &
                            avg_obs_per_file

! start of executable code

obs_num = 1
d_qc(1) = 0.0_r8

! time stored relative to jan 1st, 1900.
call set_calendar_type(GREGORIAN)
base_time = set_date(1900, 1, 1, 0, 0, 0)

!  read the necessary parameters from input.nml
call initialize_utilities()
call find_namelist_in_file("input.nml", "gtspp_to_obs_nml", iunit)
read(iunit, nml = gtspp_to_obs_nml, iostat = io)

! Record the namelist values used for the run 
if (do_output()) write(nmlfileunit, nml=gtspp_to_obs_nml)

! any needed namelist checks for sanity:

! cannot have both a single filename and a list; the namelist must
! shut one off.
if (gtspp_netcdf_file /= '' .and. gtspp_netcdf_filelist /= '') then
  call error_handler(E_ERR, 'gtspp_to_obs',                     &
                     'One of gtspp_netcdf_file or filelist must be NULL', &
                     source, revision, revdate)
endif
if (gtspp_netcdf_filelist /= '') from_list = .true.

! need to know a reasonable max number of obs that could be added here.
if (from_list) then
   call find_textfile_dims(gtspp_netcdf_filelist, nfiles, dummy)
   num_new_obs = avg_obs_per_file * nfiles
else
   num_new_obs = avg_obs_per_file
endif

! FIXME: is this a good idea?  to append to an existing file?
!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
inquire(file=gtspp_out_file, exist=file_exist)
if ( file_exist ) then

print *, "found existing obs_seq file, appending to ", trim(gtspp_out_file)
   call read_obs_seq(gtspp_out_file, 0, 0, num_new_obs, obs_seq)

else

print *, "no existing obs_seq file, creating ", trim(gtspp_out_file)
print *, "max entries = ", num_new_obs
  call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
  do k = 1, num_copies
    meta_data = 'GTSPP observation'
    call set_copy_meta_data(obs_seq, k, meta_data)
  end do
  do k = 1, num_qc
    meta_data = 'GTSPP QC'
    call set_qc_meta_data(obs_seq, k, meta_data)
  end do

end if

did_obs = .false.

! main loop that does either a single file or a list of files

filenum = 1
fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(gtspp_netcdf_filelist, filenum)
   else
      next_infile = gtspp_netcdf_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop
  
   !  open the next profile file
   call nc_check( nf90_open(next_infile, nf90_nowrite, ncid), 'file open', next_infile)

   ! time is stored in the file 2 ways: as real(double) days since 1900/1/1,
   ! and as 4 and 2 digit strings for year/mon/day/hr/min
   ! both of these are variables, not attributes

   ! start out with converting the real time.
   call nc_check( nf90_inq_varid(ncid, "time", varid) ,'inq varid time')
   call nc_check( nf90_get_var(ncid, varid, dtime)    ,'get var   time')

   ! convert to integer days and seconds, and add on to reference time.
   iday = int(dtime)
   isec = int((dtime - iday) * 86400)
   delta_time = set_time(isec, iday)
   obs_time = base_time + delta_time
   call get_time(obs_time,  osec, oday)
   
   ! get the number of depths
   call nc_check( nf90_inq_dimid(ncid, "depth", varid), 'inq dimid depth')
   call nc_check( nf90_inquire_dimension(ncid, varid, name, ndepths), 'inq dim depth')
   
   ! and read in the depth array
   call nc_check( nf90_inq_varid(ncid, "depth", varid),'inq varid depth')
   call nc_check( nf90_get_var(ncid, varid, obs_depth, &
                  start=(/1/), count=(/ndepths/)),'get var   depth')

   ! get the single lat/lon values
   call nc_check( nf90_inq_varid(ncid,"longitude",varid) ,'inq varid longitude')
   call nc_check( nf90_get_var(ncid, varid, glon),        'get var   longitude')
   call nc_check( nf90_inq_varid(ncid,"latitude",varid) ,'inq varid latitude')
   call nc_check( nf90_get_var(ncid, varid, glat),       'get var   latitude')

   ! if present, the data values from 'temperature'
   nc_rc = nf90_inq_varid(ncid,"temperature",varid) 
   if (nc_rc == nf90_noerr) then
      have_temp = .true.

      call nc_check( nf90_get_var(ncid, varid, temperature, &
            start=(/1,1,1,1/), count=(/1,1,ndepths,1/)), 'get var temperature')

      ! for now, use the data qc from netcdf file
      call nc_check( nf90_inq_varid(ncid,"TEMP_qparm",varid) ,'inq varid TEMP_qparam')
      call nc_check( nf90_get_var(ncid, varid, t_qc, &
            start=(/1,1/), count=(/1,ndepths/)),   'get var   TEMP_qparam')

   endif

   ! if present, the data values from 'salinity'
   nc_rc = nf90_inq_varid(ncid,"salinity",varid)
   if (nc_rc == nf90_noerr) then
      have_salt = .true.

      call nc_check( nf90_get_var(ncid, varid, salinity, &
            start=(/1,1,1,1/), count=(/1,1,ndepths,1/)), 'get var salinity')

      ! for now, use the data qc from netcdf file

      call nc_check( nf90_inq_varid(ncid,"PSAL_qparm",varid) ,'inq varid PSAL_qparam')
      call nc_check( nf90_get_var(ncid, varid, s_qc, &
            start=(/1,1/), count=(/1,ndepths/)),   'get var   PSAL_qparam')

   endif

   call nc_check( nf90_close(ncid) , 'close file')
   
   
 ! FIXME:
   terr = 2.0   ! temp error = 2 degrees C
   serr = 1.0   ! salinity error = 1 something?

   first_obs = .true.
   
   obsloop: do k = 1, ndepths
   
      if (have_temp) then
         ! check qc here.  if bad, skip the rest of this block
         read(t_qc(k:k), '(I1)') i_qc
         if ( i_qc /= 1 )  exit     ! in this case, 1 is good
   
         ! set qc
         d_qc(1) = 0.0    ! but for dart, a QC of 0 is good
 
         ! set location 
         if (glon < 0.0_r8) glon = glon + 180.0_r8
         call set_obs_def_location(obs_def, &
                           set_location(glon, glat, obs_depth(k),VERTISHEIGHT))
         call set_obs_def_kind(obs_def, FLOAT_TEMPERATURE)
         call set_obs_def_time(obs_def, set_time(osec, oday))
    
         call set_obs_def_error_variance(obs_def, terr * terr)
         call set_obs_def_key(obs_def, obs_num)
         call set_obs_def(obs, obs_def)
   
         obs_val(1) = temperature(k)
         call set_obs_values(obs, obs_val)
         qc_val(1)  = d_qc(1)
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

      if (have_salt) then
         ! check qc here.  if bad, skip the rest of this block
         read(s_qc(k:k), '(I1)') i_qc
         if ( i_qc /= 1 )  exit     ! in this case, 1 is good
   
         ! set qc
         d_qc(1) = 0.0    ! but for dart, a QC of 0 is good
 
         ! set location 
         if (glon < 0.0_r8) glon = glon + 180.0_r8
         call set_obs_def_location(obs_def, &
                           set_location(glon, glat, obs_depth(k),VERTISHEIGHT))
         call set_obs_def_kind(obs_def, FLOAT_SALINITY)
         call set_obs_def_time(obs_def, set_time(osec, oday))
    
         call set_obs_def_error_variance(obs_def, serr * serr)
         call set_obs_def_key(obs_def, obs_num)
         call set_obs_def(obs, obs_def)
   
         obs_val(1) = salinity(k)
         call set_obs_values(obs, obs_val)
         qc_val(1)  = d_qc(1)
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

  filenum = filenum + 1

end do fileloop

! done with main loop.  if we added any obs to the sequence, write it out.
if (did_obs) then
print *, 'ready to write, nobs = ', get_num_obs(obs_seq)
   if (get_num_obs(obs_seq) > 0) &
      call write_obs_seq(obs_seq, gtspp_out_file)

   ! minor stab at cleanup, in the off chance this will someday get turned
   ! into a subroutine in a module.  probably not all that needs to be done,
   ! but a start.
   call destroy_obs(obs)
   !call destroy_obs(prev_obs)   ! is this identical to obs?
   ! get core dumps here, not sure why?
   if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
endif

! END OF MAIN ROUTINE

! contains

! local subroutines/functions follow

! something to set the err var - make it a subroutine so we can muck with it

! subroutine to fill an obs?

end program

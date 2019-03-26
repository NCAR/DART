! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
!
! bjchoi 2014 May 15
! Because a GTSPP profile has too many data in vertical direction,
! we need to reduce the number of data in a single profile.
! Just select data near the standard (nominal) depths.

program thinned_gtspp_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     thinned_gtspp_to_obs - program that reads a list of NODS GTSPP observation 
!                    profiles of ocean temperature and salinity in netcdf
!                    format and writes a DART observation sequence file. 
!
!     created 02-sept-2009, based on the GPS reader.   nsc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use        types_mod, only : r8, MISSING_R8
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time, &
                             increment_time, get_time, set_date, operator(-),   &
                             print_date, operator(+)
use    utilities_mod, only : initialize_utilities, find_namelist_in_file,       &
                             check_namelist_read, nmlfileunit, do_output,       &
                             get_next_filename, error_handler, E_ERR, E_MSG,    &
                             find_textfile_dims, finalize_utilities
use  netcdf_utilities_mod, only : nc_check
use     location_mod, only : VERTISHEIGHT, set_location
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,         &
                             static_init_obs_sequence, init_obs, destroy_obs,   &
                             write_obs_seq, init_obs_sequence, get_num_obs,     &
                             insert_obs_in_seq, destroy_obs_sequence,           &
                             set_copy_meta_data, set_qc_meta_data, set_qc,      & 
                             set_obs_values, set_obs_def, insert_obs_in_seq
use      obs_def_mod, only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs,  &
                             set_obs_def_error_variance, set_obs_def_location,  &
                             set_obs_def_key
! FIXME: what actual instrument took these readings? FLOAT_xx is a placeholder
! for now.  must have obs_def_ocean_mod.f90 in the preprocess input list.
use     obs_kind_mod, only : QTY_TEMPERATURE, QTY_SALINITY,                   &
                             FLOAT_TEMPERATURE, FLOAT_SALINITY

use           netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=256) :: next_infile
character (len=NF90_MAX_NAME) :: dimname
integer :: ncid, varid, ndepths, k, nfiles, num_new_obs,  &
           oday, osec,                   &
           iday, isec,               &
           obs_num, io, iunit, filenum, dummy, i_qc, nc_rc
logical :: file_exist, first_obs, did_obs, from_list = .false.
logical :: have_temp, have_salt
real(r8) :: terr, serr,         & 
            obs_val(1), qc_val(1),  &
            dtime, glat, glon, d_qc(1)

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, base_time, delta_time

! initialize some values
integer, parameter :: nmaxdepths = 5000   !  max number of observation depths
real(r8) :: obs_depth(nmaxdepths)   = -1.0_r8
real(r8) :: temperature(nmaxdepths) = MISSING_R8
real(r8) :: salinity(nmaxdepths)    = MISSING_R8
integer :: t_qc(nmaxdepths), s_qc(nmaxdepths)

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=128) :: gtspp_netcdf_file     = '1234567.nc'
character(len=128) :: gtspp_netcdf_filelist = 'gtspp_to_obs_filelist'
character(len=128) :: gtspp_out_file        = 'obs_seq.gtspp'
integer            :: avg_obs_per_file      = 500
logical            :: debug                 = .false.

namelist /thinned_gtspp_to_obs_nml/ gtspp_netcdf_file, &
                                    gtspp_netcdf_filelist, &
                                    gtspp_out_file, &
                                    avg_obs_per_file, &
                                    debug

!pseudo standard depth for a profile data. 5 lines were added by bjchoi
integer  :: mdepth
real(r8) :: std_depth(22)
DATA std_depth(1:22) / 0.0, 15.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, & 
                     400.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, & 
                    4000.0, 5000.0, 6000.0, 8000.0, 12000.0 /


! start of executable code

obs_num = 1
d_qc(1) = 0.0_r8

! time stored relative to jan 1st, 1900.
call set_calendar_type(GREGORIAN)
base_time = set_date(1900, 1, 1, 0, 0, 0)

!  read the necessary parameters from input.nml
call initialize_utilities()
call find_namelist_in_file("input.nml", "thinned_gtspp_to_obs_nml", iunit)
read(iunit, nml = thinned_gtspp_to_obs_nml, iostat = io)

! Record the namelist values used for the run 
if (do_output()) write(nmlfileunit, nml=thinned_gtspp_to_obs_nml)

! any needed namelist checks for sanity:

! cannot have both a single filename and a list; the namelist must
! shut one off.
if (gtspp_netcdf_file /= '' .and. gtspp_netcdf_filelist /= '') then
  call error_handler(E_ERR, 'thinned_gtspp_to_obs',                     &
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
! for now, disable this.  too easy to make files with duplicate obs in them.
if ( file_exist ) then
!if (.false.) then

print *, "found existing obs_seq file, appending to ", trim(gtspp_out_file)
   call read_obs_seq(gtspp_out_file, 0, 0, num_new_obs, obs_seq)

else

print *, "no existing obs_seq file, creating ", trim(gtspp_out_file)
print *, "max entries = ", num_new_obs
  call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
  do k = 1, num_copies
    call set_copy_meta_data(obs_seq, k, 'GTSPP observation')
  end do
  do k = 1, num_qc
    call set_qc_meta_data(obs_seq, k, 'GTSPP QC')
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
  
   have_temp = .false.
   have_salt = .false.

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
   call nc_check( nf90_inq_dimid(ncid, "z", varid), 'inq dimid depth')
   call nc_check( nf90_inquire_dimension(ncid, varid, len=ndepths, name=dimname), 'inq dim depth')
   
   ! and read in the depth array
   call nc_check( nf90_inq_varid(ncid, "z", varid),'inq varid depth')
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

      ! for now, use the data qc from netcdf file 't_qc'
      call nc_check( nf90_inq_varid(ncid,"temperature_quality_flag",varid) ,'inq varid TEMP_qparam')
      call nc_check( nf90_get_var(ncid, varid, t_qc, &
            start=(/1/), count=(/ndepths/)),   'get var   TEMP_qparam')

   endif

   ! if present, the data values from 'salinity'
   nc_rc = nf90_inq_varid(ncid,"salinity",varid)
   if (nc_rc == nf90_noerr) then
      have_salt = .true.

      call nc_check( nf90_get_var(ncid, varid, salinity, &
            start=(/1,1,1,1/), count=(/1,1,ndepths,1/)), 'get var salinity')

      ! for now, use the data qc from netcdf file 's_qc'

      call nc_check( nf90_inq_varid(ncid,"salinity_quality_flag",varid) ,'inq varid PSAL_qparam')
      call nc_check( nf90_get_var(ncid, varid, s_qc, &
            start=(/1/), count=(/ndepths/)),   'get var   PSAL_qparam')

   endif

   call nc_check( nf90_close(ncid) , 'close file')
   
   if (debug) then
      print *, 'input file: ', trim(next_infile)
      print *, 'ndepths: ', ndepths
      print *, 'has temps, salinity: ', have_temp, have_salt
      call print_date(obs_time, 'obs time')
   endif
   
 ! FIXME: these have no physical basis; selected to get us running
 ! but need some domain expertise for guidance.
   terr = 0.5     ! temp error = fixed at 1 degrees C
   serr = 0.05    ! salinity error = 1 g/kg, which is 0.001 kg/kg

   first_obs = .true.
   mdepth = 1
   
   obsloop: do k = 1, ndepths
   
    ! 5 lines were ADDED by bjchoi 2015.05.15.
    ! It could be improved better than this.
    if( obs_depth(k) .GE. std_depth(mdepth) ) then
       mdepth = mdepth + 1 
    else
       t_qc(k) = -1
       s_qc(k) = -1
    endif    

     ! check qc here.  if bad, skip the rest of this block
     i_qc = t_qc(k)
     if (have_temp .and. i_qc == 1) then

         ! set qc to a good dart val
         d_qc(1) = 0.0    ! for dart, a QC of 0 is good
 
         ! set location - incoming obs are -180 to 180 in longitude;
         ! dart wants 0 to 360.
         if (glon < 0.0_r8) glon = glon + 360.0_r8
         call set_obs_def_location(obs_def, &
                           set_location(glon, glat, obs_depth(k),VERTISHEIGHT))
         call set_obs_def_type_of_obs(obs_def, FLOAT_TEMPERATURE)
         call set_obs_def_time(obs_def, set_time(osec, oday))
    
         call set_obs_def_error_variance(obs_def, terr * terr)
         call set_obs_def_key(obs_def, obs_num)
         call set_obs_def(obs, obs_def)
   
         obs_val(1) = temperature(k)
         call set_obs_values(obs, obs_val)
         qc_val(1)  = d_qc(1)
         call set_qc(obs, qc_val)
    
         ! first one, insert with no prev.  otherwise, since all times are the
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


      ! check qc here.  if bad, skip the rest of this block
      i_qc = s_qc(k)
      if (have_salt .and. i_qc == 1) then

         ! set qc to good dart val
         d_qc(1) = 0.0    ! but for dart, a QC of 0 is good
 
         ! set location - incoming obs are -180 to 180 in longitude;
         ! dart wants 0 to 360.
         if (glon < 0.0_r8) glon = glon + 360.0_r8
         call set_obs_def_location(obs_def, &
                           set_location(glon, glat, obs_depth(k),VERTISHEIGHT))
         call set_obs_def_type_of_obs(obs_def, FLOAT_SALINITY)
         call set_obs_def_time(obs_def, set_time(osec, oday))
    
         call set_obs_def_error_variance(obs_def, serr * serr)
         call set_obs_def_key(obs_def, obs_num)
         call set_obs_def(obs, obs_def)
   
         ! incoming obs are g/kg (practical salinity units - psu)
         ! model works in kg/kg (model salinity units - msu) so convert here.
         ! obs_val(1) = salinity(k) / 1000.0_r8
         obs_val(1) = salinity(k)  ! ROMS use psu            
         call set_obs_values(obs, obs_val)
         qc_val(1)  = d_qc(1)
         call set_qc(obs, qc_val)
    
         ! first one, insert with no prev.  otherwise, since all times are the
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

call error_handler(E_MSG,'thinned_gtspp_to_obs','Finished successfully.',source,revision,revdate)
call finalize_utilities()

! END OF MAIN ROUTINE

! contains

! local subroutines/functions follow

! something to set the err var - make it a subroutine so we can muck with it

! subroutine to fill an obs?

end program thinned_gtspp_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

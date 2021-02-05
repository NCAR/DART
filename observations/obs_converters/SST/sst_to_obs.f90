! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program sst_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     sst_to_obs - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, MISSING_R8

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time, &
                              increment_time, get_time, set_date, operator(-),   &
                              print_date, operator(+)

use     utilities_mod, only : initialize_utilities, find_namelist_in_file,       &
                              check_namelist_read, nmlfileunit, do_output,       &
                              get_next_filename, error_handler, E_ERR, E_MSG,    &
                              find_textfile_dims, finalize_utilities,            &
                              do_nml_file, do_nml_term

use  netcdf_utilities_mod, only : nc_check

use      location_mod, only : VERTISSURFACE, set_location

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,         &
                              static_init_obs_sequence, init_obs, destroy_obs,   &
                              write_obs_seq, init_obs_sequence, get_num_obs,     &
                              insert_obs_in_seq, destroy_obs_sequence,           &
                              set_copy_meta_data, set_qc_meta_data, set_qc,      & 
                              set_obs_values, set_obs_def, insert_obs_in_seq

use       obs_def_mod, only : obs_def_type, &
                              set_obs_def_key, &
                              set_obs_def_time, &
                              set_obs_def_type_of_obs, &
                              set_obs_def_error_variance, &
                              set_obs_def_location

use      obs_kind_mod, only : QTY_TEMPERATURE, SATELLITE_INFRARED_SST 

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   '$URL$'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

character(len=256) :: next_infile
character(len=512) :: string1, string2

integer :: ncid, varid, io, iunit, filenum
integer :: oday, osec, iday, isec
integer :: num_new_obs, obs_num
integer :: i, j, k, nlat, nlon
integer :: i_base, j_base, tmp_day

logical :: first_obs, from_list = .false.

real(r8) :: obs_val(1), qc_val(1), d_qc(1), dtime

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, base_time, delta_time

real(r8), allocatable :: temperature(:,:)
real(r8), allocatable :: t_err(:,:)
real(r8), allocatable :: glat(:), glon(:)
real(r8) :: scale_factor, add_offset
real(r8) :: tmp_qc

integer, allocatable  ::   t_qc(:,:)
integer, allocatable  :: t_mask(:,:)

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=256) :: sst_netcdf_file     = '1234567.nc'
character(len=256) :: sst_netcdf_filelist = 'sst_to_obs_filelist'
character(len=256) :: sst_out_file        = 'obs_seq.sst'
logical            :: debug               = .false.
integer            :: subsample_intv      = 1

! representation error, which depends on model/data grid size, 
! is a very small value for NOAA OI SST
real(r8) :: sst_rep_error = 0.3  ! minimum SST observation error

namelist /sst_to_obs_nml/ sst_netcdf_file,                   &
                          sst_netcdf_filelist, sst_out_file, &
                          debug, subsample_intv, sst_rep_error

! start of executable code


! time is stored relative to Jan 1, 1981 for NOAA OI sst data.

call set_calendar_type(GREGORIAN)
base_time = set_date(1981, 1, 1, 0, 0, 0)

!  read the necessary parameters from input.nml
call initialize_utilities()
call find_namelist_in_file('input.nml', 'sst_to_obs_nml', iunit)
read(iunit, nml = sst_to_obs_nml, iostat = io)

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=sst_to_obs_nml)
if (do_nml_term()) write(     *     , nml=sst_to_obs_nml)

! cannot have both a single filename and a list;
! the namelist must shut one off.

if (sst_netcdf_file /= '' .and. sst_netcdf_filelist /= '') then
  call error_handler(E_ERR, 'sst_to_obs',                     &
                     'One of sst_netcdf_file or filelist must be NULL', &
                     source, revision, revdate)
endif
if (sst_netcdf_filelist /= '') from_list = .true.

! Get number of observations
if (from_list) then
   next_infile = get_next_filename(sst_netcdf_filelist, 1)
else
   next_infile = sst_netcdf_file
endif

call nc_check( nf90_open(next_infile, nf90_nowrite, ncid), 'open '//trim(next_infile))
call nc_check( nf90_inq_dimid(ncid, 'lon', varid), 'inq dimid lon')
call nc_check( nf90_inquire_dimension(ncid, varid, len=nlon), 'inq dimlon')
call nc_check( nf90_inq_dimid(ncid, 'lat', varid), 'inq dimid lat')
call nc_check( nf90_inquire_dimension(ncid, varid, len=nlat), 'inq dimlat')
call nc_check( nf90_close(ncid) , 'close file')

!>@todo FIXME ... num_new_obs is only correct if there is only 1 input file
num_new_obs = nlon*nlat

! Initialize
allocate(temperature(nlon,nlat))
allocate(      t_err(nlon,nlat))
allocate(       t_qc(nlon,nlat))
allocate(     t_mask(nlon,nlat))

allocate(glat(nlat), glon(nlon))

call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)

do k = 1, num_copies
  call set_copy_meta_data(obs_seq, k, 'SST observation')
enddo

do k = 1, num_qc
  call set_qc_meta_data(obs_seq, k, 'SST QC')
enddo

obs_num = 1
d_qc(1) = 0.0_r8
first_obs = .true.
   
! main loop that does either a single file or a list of files

filenum = 1
fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(sst_netcdf_filelist, filenum)
   else
      next_infile = sst_netcdf_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop
  
   !  open the next profile file
   call nc_check( nf90_open(next_infile, nf90_nowrite, ncid), 'file open', next_infile)

   ! time is stored in the file 2 ways: as real(double) seconds since 1981/1/1,
   ! and as 4 and 2 digit strings for year/mon/day/hr/min
   ! both of these are variables, not attributes

   ! start out with converting the real time.
   call nc_check( nf90_inq_varid(ncid, 'time', varid) ,'inq varid time')
   call nc_check( nf90_get_var(ncid, varid, dtime)    ,'get var   time')

   ! convert to integer days and seconds, and add on to reference time.
   iday = int( dtime / 86400.0_r8)
   isec = int( dtime - iday * 86400)
   delta_time = set_time(isec, iday)
   obs_time = base_time + delta_time
   call get_time(obs_time,  osec, oday)
   
   ! get the lat/lon arrays

   io = nf90_inq_varid(ncid, 'lon', varid)
   call nc_check(io,'inq_varid "lon"')

   io = nf90_get_var(ncid, varid, glon)
   call nc_check(io, 'get_var "lon"')

   ! ensure longitudes are [0,360] 
   where(glon < 0.0_r8) glon = glon + 360.0_r8

   io = nf90_inq_varid(ncid, 'lat', varid)
   call nc_check(io, 'inq_varid "lat"')

   io = nf90_get_var(ncid, varid, glat)
   call nc_check(io, 'get_var "lat"')
  
   ! if present, the data values from 'temperature'
   io = nf90_inq_varid(ncid, 'analysed_sst', varid) 

   if (io /= nf90_noerr ) then
      call error_handler(E_MSG, 'sst_to_obs', trim(next_infile)//' has no sst') 
      cycle fileloop
   endif

   ! Get temperature, scale and offset

   io = nf90_get_var(ncid, varid, temperature, &
                     start=(/1,1,1/), count=(/nlon,nlat,1/))
   call nc_check(io, 'get_var analysed_sst')

   io = nf90_get_att(ncid, varid, 'scale_factor',scale_factor)
   call nc_check(io, 'get_att analysed_sst scale_factor')

   io = nf90_get_att(ncid, varid, 'add_offset', add_offset)
   call nc_check(io, 'get_att analysed_sst add_offset')
  
   temperature = temperature * scale_factor + add_offset

   ! Get analysis_error, scale and offset
 
   io = nf90_inq_varid(ncid,'analysis_error',varid)
   call nc_check(io,'inq_varid analysis_error')

   io = nf90_get_var(ncid, varid, t_qc, &
                     start=(/1,1,1/), count=(/nlon,nlat,1/))
   call nc_check(io, 'get_var analysis_error')

   io = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
   call nc_check(io, 'get att analysis_error scale_factor')

   io = nf90_get_att(ncid, varid, 'add_offset', add_offset)
   call nc_check(io, 'get att analysis_error add_offset')

   t_qc = t_qc * scale_factor + add_offset

   ! Get data mask

   io = nf90_inq_varid(ncid, 'mask', varid)
   call nc_check(io, 'inq_varid mask')

   io = nf90_get_var(ncid, varid, t_mask, start=(/1,1,1/), count=(/nlon,nlat,1/))
   call nc_check(io, 'get_var mask')

   call nc_check( nf90_close(ncid) , 'closing "'//trim(next_infile)//'"')
   
   tmp_day = mod( oday,subsample_intv*subsample_intv )
   if( tmp_day == 0 ) tmp_day = subsample_intv*subsample_intv
   j_base = 1 + int((tmp_day-1)/subsample_intv)
   i_base = mod(tmp_day,subsample_intv) 
   if( i_base == 0 ) i_base = subsample_intv

   obslooplat: do j = j_base, nlat, subsample_intv
   obslooplon: do i = i_base, nlon, subsample_intv
   
     ! check qc here.  if bad, skip the rest of this block
     ! we want small observation error smaller than 1 degree C
     tmp_qc = t_qc(i,j) + sst_rep_error

     if (tmp_qc < 20.0_r8 .and. t_mask(i,j) == 1.0_r8) then

         ! set qc to a good dart val
         d_qc(1) = 0.0    ! for dart, a QC of 0 is good

         call set_obs_def_location(obs_def, &
                         set_location(glon(i), glat(j), 0.0_r8, VERTISSURFACE))
         call set_obs_def_type_of_obs(obs_def, SATELLITE_INFRARED_SST)
         call set_obs_def_time(obs_def, set_time(osec, oday))
         call set_obs_def_error_variance(obs_def, tmp_qc * tmp_qc)
         call set_obs_def_key(obs_def, obs_num)
         call set_obs_def(obs, obs_def)
   
         obs_val(1) = temperature(i,j) - 273.15_r8
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

         if (obs_num > num_new_obs) then
            write(string1,*)'reached maximum number of observations while'
            write(string2,*)'reading file # ',filenum,' "'//trim(next_infile)//'"'
            call error_handler(E_ERR, 'sst_to_obs', string1, &
                       source, revision, revdate, text2=string2)
         endif
 
      endif

   enddo obslooplon
   enddo obslooplat

  filenum = filenum + 1

end do fileloop

! done with main loop.  if we added any obs to the sequence, write it out.
if (obs_num > 0) then
   print *, 'ready to write, nobs = ', get_num_obs(obs_seq)
   if (get_num_obs(obs_seq) > 0) &
      call write_obs_seq(obs_seq, sst_out_file)

   ! minor stab at cleanup, in the off chance this will someday get turned
   ! into a subroutine in a module.  probably not all that needs to be done,
   ! but a start.
   call destroy_obs(obs)
   !call destroy_obs(prev_obs)   ! is this identical to obs?
   ! get core dumps here, not sure why?
   if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
else
   call error_handler(E_ERR,'sst_to_obs','no obs converted', &
              source, revision,revdate)
endif

call error_handler(E_MSG,'sst_to_obs','Finished successfully.')
call finalize_utilities()

end program sst_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

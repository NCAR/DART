! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program netCDF_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! netCDF_to_obs
!
! This converter is intended to be used as example of how to convert
! basic netCDF files to observation sequence files in one step.
! The netCDF file(s) in question is/are (bi)monthly values, 
! so we want to put each time in a separate observation sequence file.
! The header of the (hypothetical) netCDF file looks like this:
!
! ncdump -h GIMMS3g_LAI_Bimonthly_2000_2010_60min.nc
!     netcdf GIMMS3g_LAI_Bimonthly_2000_2010_60min {
!     dimensions:
!         time = 264 ;
!         lat = 180 ;
!         lon = 360 ;
!         nv = 2 ;
!     variables:
!         float LAI(time, lat, lon) ;
!            LAI:_FillValue = -9999.f ;
!            LAI:units = "m^2/m^2" ;
!            LAI:cell_methods = "time: mean" ;
!            LAI:standard_name = "leaf_area_index" ;
!            LAI:long_name = "Leaf Area Index" ;
!         double lat(lat) ;
!            lat:standard_name = "latitude" ;
!            lat:long_name = "latitude coordinate" ;
!            lat:units = "degrees_north" ;
!         double lon(lon) ;
!            lon:standard_name = "longitude" ;
!            lon:long_name = "longitude coordinate" ;
!            lon:units = "degrees_east" ;
!         double time(time) ;
!            time:units = "days since 2000-01-01 00:00:00" ;
!            time:calendar = "gregorian" ;
!            time:standard_name = "time" ;
!            time:climatology = "climatology_bounds" ;
!         int climatology_bounds(time, nv) ;
!            climatology_bounds:calendar = "standard" ;
!            climatology_bounds:description = "start and end time for each time stamp" ;
!            climatology_bounds:units = "days since 2000-01-01 00:00:00" ;
!     
!     // global attributes:
!         :source = "AVHRR GIMMS LAI3g version 2" ;
!}
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : digits12, r8, MISSING_R8

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_time, get_time, set_date, get_date,  &
                              print_date, print_time

use     utilities_mod, only : initialize_utilities, finalize_utilities,   &
                              check_namelist_read, find_namelist_in_file, &
                              do_nml_file, nmlfileunit,                   &
                              do_nml_term, logfileunit,                   &
                              do_output, get_next_filename,               &
                              error_handler, E_ERR, E_MSG, E_WARN

use netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, &
                                 nc_get_dimension_size,                &
                                 nc_variable_exists,                   &
                                 nc_get_variable_info,                 &
                                 nc_get_variable,                      &
                                 nc_get_attribute_from_variable

use  obs_utilities_mod, only : getvar_real_2d, create_3d_obs

use       location_mod, only : VERTISUNDEF, set_location

use  obs_sequence_mod, only : obs_type, obs_sequence_type,             &
                              static_init_obs_sequence,                &
                              init_obs, destroy_obs,                   &
                              init_obs_sequence, destroy_obs_sequence, &
                              set_copy_meta_data, set_qc_meta_data,    &
                              insert_obs_in_seq,                       &
                              get_num_obs, write_obs_seq

use      obs_kind_mod, only : get_index_for_type_of_obs

use netcdf

implicit none

character(len=*), parameter :: source   = 'netCDF_to_obs.f90'

integer, parameter :: NUM_COPIES = 1
integer, parameter :: NUM_QC     = 1

character(len=256) :: next_infile
character(len=512) :: string1, string2, string3

integer :: ncid, io, iunit, filenum
integer :: day, seconds
integer :: num_new_obs
integer :: i, j, nlat, nlon, itime
integer :: obs_type_integer

logical :: first_obs, from_list = .false.

real(r8) :: qc, lat, lon, vertvalue
real(r8) :: missing_value

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time

real(r8), allocatable :: Latitude(:), Longitude(:)
real(r8), allocatable :: observation(:,:)
real(digits12), allocatable :: obstime(:)

integer :: ndims, idimension, ntimes, timedimid
integer                      :: ncstart( NF90_MAX_VAR_DIMS)
integer                      :: nccount( NF90_MAX_VAR_DIMS)
integer                      :: dimlens( NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: dimnames(NF90_MAX_VAR_DIMS)
character(len=256) :: output_file

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=256) :: input_file      = ''
character(len=256) :: input_file_list = ''
character(len=256) :: output_file_base = 'obs_seq.out'
character(len=256) :: observation_varname
character(len=256) :: observation_type
logical            :: debug           = .false.
real(r8)           :: lon_bounds(2) = (/   0.0_r8, 360.0_r8 /)
real(r8)           :: lat_bounds(2) = (/ -90.0_r8,  90.0_r8 /)

! must include representation error, which depends on model/data grid size,
real(r8) :: obs_error_standard_deviation = 0.2

namelist /netCDF_to_obs_nml/ &
   input_file, input_file_list, output_file_base, &
   observation_varname, observation_type, &
   debug, obs_error_standard_deviation, &
   lon_bounds, lat_bounds

call set_calendar_type(GREGORIAN)

!  read the necessary parameters from input.nml
call initialize_utilities(source)
call static_init_obs_sequence()

call find_namelist_in_file('input.nml', 'netCDF_to_obs_nml', iunit)
read(iunit, nml = netCDF_to_obs_nml, iostat = io)
call check_namelist_read(iunit, io, 'netCDF_to_obs_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=netCDF_to_obs_nml)
if (do_nml_term()) write(     *     , nml=netCDF_to_obs_nml)

! cannot have both a single filename and a list;
! the namelist must shut one off.

if (input_file /= '' .and. input_file_list /= '') then
   call error_handler(E_ERR, 'main', 'One of input_file or filelist must be NULL', &
                     source)
endif

if (input_file_list /= '') from_list = .true.

!-----------------------------------------------------------------------
! Match the string defining the observation_type to the integer code

obs_type_integer = get_index_for_type_of_obs(observation_type)
if (obs_type_integer < 0) then
   write(string1,*)'observation_type "'//trim(observation_type)//'" is not a valid TYPE.'
   write(string2,*)'valid TYPEs are summarized in obs_kind_mod.f90 after processing by preprocess.'
   write(string3,*)'They are defined in the obs_def_XXXXX_mod.f90 files.'
   call error_handler(E_ERR, 'main', string1, source, text2=string2, text3=string3)
endif

!-----------------------------------------------------------------------
! main loop that does either a single file or a list of files
! Each timestep is going to a separate output file.

filenum = 1
fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then
      next_infile = get_next_filename(input_file_list, filenum)
   else
      next_infile = input_file
      if (filenum > 1) next_infile = ''
   endif

   write(*,*)'filenum ',filenum,', is next_infile ="'//trim(next_infile)//'"'
   if (next_infile == '') exit fileloop

   ncid  = nc_open_file_readonly(next_infile, 'netCDF_to_obs')

   call nc_get_variable_info(ncid, observation_varname, &
                             ndims   = ndims,   &
                             dimlens = dimlens, &
                             dimnames= dimnames, &
                             context = 'netCDF_to_obs setup')

   ! Each timestep will be written to a separate output file whose
   ! name reflects the time. All other dimensions are used to calculate
   ! the theoretical maximum number of observations in each file.

   num_new_obs = 1
   do idimension = 1,ndims
      select case (dimnames(idimension))
         case ('time','Time','TIME')
            ntimes = dimlens(idimension)
            ncstart(idimension) = 1
            nccount(idimension) = 1
            timedimid = idimension
         case default
            num_new_obs = num_new_obs * dimlens(idimension)
            ncstart(idimension) = 1
            nccount(idimension) = dimlens(idimension)
            continue
      end select
   enddo

   allocate(obstime(ntimes))

   call nc_read_time(ncid, obstime)

   nlat  = nc_get_dimension_size(ncid,'lat')
   nlon  = nc_get_dimension_size(ncid,'lon')

   allocate(Latitude(nlat))
   allocate(Longitude(nlon))
   allocate(observation(nlon,nlat))  ! just for one slab

   call nc_get_variable(ncid, 'lon', Longitude)
   call nc_get_variable(ncid, 'lat',  Latitude)

   ! ensure longitudes are [0,360]
   where(Longitude < 0.0_r8) Longitude = Longitude + 360.0_r8

   timeloop : do itime = 1,ntimes

      if (debug) write(*,*)'timestep ',itime,' of ',ntimes

      first_obs = .true.

      call init_obs(              obs, NUM_COPIES, NUM_QC)
      call init_obs(         prev_obs, NUM_COPIES, NUM_QC)
      call init_obs_sequence( obs_seq, NUM_COPIES, NUM_QC, num_new_obs)
      call set_copy_meta_data(obs_seq, NUM_COPIES, 'observation')
      call set_qc_meta_data(  obs_seq, NUM_QC,     'original QC')

      day      = floor(obstime(itime))
      seconds  = nint((obstime(itime) - real(day,digits12)) * 86400.0_digits12)
      obs_time = set_time(seconds,day)

      ncstart(timedimid) = itime

      if (debug) then
         write(*,*)'ncstart is ',ncstart(1:ndims)
         write(*,*)'nccount is ',nccount(1:ndims)
      endif

      call getvar_real_2d(ncid, observation_varname, observation, dmiss=missing_value, &
                           nc_start=ncstart(1:ndims), nc_count=nccount(1:ndims))

      latloop: do j = 1, nlat
         lat = Latitude(j)
         if ( lat < lat_bounds(1) ) cycle latloop
         if ( lat > lat_bounds(2) ) cycle latloop
         lonloop: do i = 1, nlon

            lon = Longitude(i)
            if ( lon < lon_bounds(1) ) cycle lonloop
            if ( lon > lon_bounds(2) ) cycle lonloop

            if (observation(i,j) == missing_value ) cycle lonloop

            !------------------------------------------------------------------
            ! Create observation
            

            vertvalue  = -888888.0_r8  ! does not matter for LAI
            qc         = 0       ! they are all good

            call create_3d_obs(lat, lon, vertvalue, VERTISUNDEF, &
                         observation(i,j), obs_type_integer, &
                         obs_error_standard_deviation, day, seconds, qc, obs)

            ! The insert operation is fast if given a time-ordered starting point.
            ! Since all these observations are at the same time, just use the
            ! previous observation as the time ordered starting point.
            if (first_obs) then
               call insert_obs_in_seq(obs_seq, obs)
               first_obs = .false.
            else
              call insert_obs_in_seq(obs_seq, obs, prev_obs)
            endif
            prev_obs = obs

         enddo lonloop
      enddo latloop

      ! if we added any obs to the sequence, write it out.
      if (get_num_obs(obs_seq) > 0) then

         call construct_filename(output_file_base, obs_time, output_file)

         write(string1,*) 'writing ', get_num_obs(obs_seq), ' observations to "'//trim(output_file)//'"'
         call error_handler(E_MSG, 'main', string1, source)

         call write_obs_seq(obs_seq, output_file)

      else
         write(string1,*) ' no observations for timestep ', itime
         write(string2,*) ' in file "'//trim(next_infile)//'"'
         call print_date(obs_time,' no observations at this date')
         call print_time(obs_time,' no observations at this time')
         call error_handler(E_WARN, 'main', string1, source, text2=string2)
      endif

      ! Little known fact that if obs = prev_obs, destroying one of them destroys them both
      call destroy_obs(obs)
      ! if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
      call destroy_obs_sequence(obs_seq)

      ! a little whitespace helps output readability
      write(*,*)''
      write(logfileunit,*)''

   end do timeloop

   deallocate(Latitude)
   deallocate(Longitude)
   deallocate(observation)

   call nc_close_file(ncid)
   filenum = filenum + 1

enddo fileloop


call error_handler(E_MSG, 'main','Finished successfully.',source)
call finalize_utilities()

contains

!-----------------------------------------------------------------------
!> read the time variable from a netCDF file, applying 'common' attributes
!> like the 'days since ....' and 'calendar' ...
!> If the variable is not named 'time', 'Time', or 'TIME', the optional
!> argument may be used to specify the variable name.
!>

subroutine nc_read_time(ncid, timearray, varname)

integer,                    intent(in)  :: ncid
real(digits12),             intent(out) :: timearray(:)
character(len=*), optional, intent(in)  :: varname

character(len=NF90_MAX_NAME) :: timevariablename
character(len=256)           :: unitstring
type(time_type)              :: calendar_start
integer                      :: ios, year, month, days, hours, minutes, seconds
real(digits12)               :: timeorigin
integer                      :: ntimes

! Search for the time variable name.
write(string1,*)'time variable is not one of "time", "Time", or "TIME"'
if (present(varname)) then
   timevariablename = adjustl(varname)
   if (.not. nc_variable_exists(ncid,timevariablename)) then
      write(string1,*)'time variable is not "'//trim(timevariablename)//'"'
      call error_handler(E_ERR,'nc_read_time',string1,source,text2='variable not found')
   endif
elseif (nc_variable_exists(ncid,'time')) then
   timevariablename =         'time'
elseif (nc_variable_exists(ncid,'Time')) then
   timevariablename =         'Time'
elseif (nc_variable_exists(ncid,'TIME')) then
   timevariablename =         'TIME'
else
   call error_handler(E_ERR,'nc_read_time',string1,source,text2='variable not found')
endif

! FIXME set the calendar, check for legitimate values
! calendar = nc_get_calendar(ncid,timevariablename)
!
! call error_handler(E_ERR, 'nc_read_time', &
!     'calendar type "'//trim(calendar)//'" unsupported by nc_read_time() routine', &
!                       source, text2=string3)

call nc_get_variable(ncid,timevariablename,timearray,'nc_read_time')
call nc_get_attribute_from_variable(ncid, timevariablename, 'units', unitstring, 'nc_read_time')

if (unitstring(1:10) == 'days since') then

   read(unitstring,'(11x,i4,5(1x,i2))',iostat=ios)year,month,days,hours,minutes,seconds
   if (ios /= 0) then
      write(string1,*)'Unable to interpret time "units" attribute. Error status was ',ios
      write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS", got "'//trim(unitstring)//'"'
      call error_handler(E_ERR, 'nc_read_time:', string1, source, text2=string2)
   endif

else if (unitstring(1:13) == 'seconds since') then

   read(unitstring,'(14x,i4,5(1x,i2))',iostat=ios)year,month,days,hours,minutes,seconds
   if (ios /= 0) then
      write(string1,*)'Unable to interpret time "units" attribute. Error status was ',ios
      write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS", got "'//trim(unitstring)//'"'
      call error_handler(E_ERR, 'read_model_time:', string1, source, text2=string2)
   endif

else

   write(string1,*) 'looking for "days since" or "seconds since" in the "units" attribute'
   write(string2,*) 'units read as "'//trim(unitstring)//'"'
   call error_handler(E_ERR, 'read_model_time:', 'unable to set start of calendar', &
                      source, text2=string1, text3=string2)

endif

! set_date is relative to the calendar in use already
calendar_start = set_date(year, month, days, hours, minutes, seconds)

call get_time(calendar_start, seconds, days)

timeorigin = days + seconds/86400.0_r8
timearray  = timearray + timeorigin

! Just confirm that the first time is what we expect.
! ignore the reuse of 'calendar_start'
if ( debug ) then

   days       = floor(timearray(1))
   seconds    = nint((timearray(1) - real(days,digits12)) * 86400.0_digits12)
   calendar_start = set_time(seconds,days)

   call print_time(calendar_start,'First time in file is')
   call print_date(calendar_start,'First date in file is')

   ntimes = size(timearray)

   days       = floor(timearray(ntimes))
   seconds    = nint((timearray(ntimes) - real(days,digits12)) * 86400.0_digits12)
   calendar_start = set_time(seconds,days)

   call print_time(calendar_start,'Last  time in file is')
   call print_date(calendar_start,'Last  date in file is')

endif

end subroutine nc_read_time


!-----------------------------------------------------------------------
!> Append the CESM-style date (YYYY-MM-DD-sssss) to the filename base


subroutine construct_filename(base, mytime, filename)

character(len=*), intent(in)  :: base
type(time_type),  intent(in)  :: mytime
character(len=*), intent(out) :: filename

integer :: year, month, days, hours, minutes, seconds

call get_date(mytime, year, month, days, hours, minutes, seconds)

seconds = (hours*60 + minutes)*60 + seconds

write(filename,100) trim(base),year,month,days,seconds

100 format(A,'.',i4.4,'-',i2.2,'-',i2.2,'-',i5.5)

end subroutine construct_filename


end program netCDF_to_obs

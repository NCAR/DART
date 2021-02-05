! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program oi_sst_to_obs

!> title = "NOAA High-resolution Blended Analysis: Daily Values using AVHRR only" ;
!> institution = "NOAA/NCDC" ;
!> source = "NOAA/NCDC  ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/" ;
!> comment = "Reynolds, et al., 2007: 
!>            Daily High-Resolution-Blended Analyses for Sea Surface Temperature.
!>            J. Climate, 20, 5473-5496.
!>            Climatology is based on 1971-2000 OI.v2 SST,
!>            Satellite data: Navy NOAA17 NOAA18 AVHRR, Ice data: NCEP ice." ;
!> references = "https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html" ;
!> dataset_title = "NOAA Daily Optimum Interpolation Sea Surface Temperature" ;

! 	float sst(time, lat, lon) ;
!            sst:long_name = "Daily Sea Surface Temperature" ;
!            sst:units = "degC" ;
!            sst:valid_range = -3.f, 45.f ;
!            sst:missing_value = -9.96921e+36f ;
!            sst:precision = 2.f ;
!            sst:dataset = "NOAA High-resolution Blended Analysis" ;
!            sst:var_desc = "Sea Surface Temperature" ;
!            sst:level_desc = "Surface" ;
!            sst:statistic = "Mean" ;
!            sst:parent_stat = "Individual Observations" ;
!            sst:actual_range = -1.8f, 35.54f ;

use         types_mod, only : r8, digits12

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_time, get_time, print_time, &
                              set_date, get_date, print_date, &
                              operator(+), operator(-)

use     utilities_mod, only : initialize_utilities, find_namelist_in_file, &
                              check_namelist_read, nmlfileunit, &
                              error_handler, E_ERR, E_MSG, &
                              finalize_utilities, do_nml_file, do_nml_term

use      location_mod, only : VERTISSURFACE, set_location

use  obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs, &
                              static_init_obs_sequence, init_obs_sequence, &
                              set_copy_meta_data, set_qc_meta_data, &
                              get_num_obs, write_obs_seq, destroy_obs_sequence

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use      obs_kind_mod, only : OI_SEA_SURFACE_TEMPERATURE 

use netcdf_utilities_mod, only : nc_check, nc_open_file_readonly, nc_close_file, &
                                 nc_get_variable, nc_get_attribute_from_variable, &
                                 nc_get_dimension_size, nc_get_variable_size 

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   '$URL$'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'
character(len=*), parameter :: routine  = 'oi_sst_to_obs'

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

real(r8), parameter :: qc = 0.0_r8

character(len=256) :: output_file
character(len=512) :: string1, string2

integer :: ncid, varid, io, iunit
integer :: oday, osec, iday, isec
integer :: year, month, day, hour, minutes, seconds
integer :: num_new_obs, nmissing
integer :: i, j, nlat, nlon, ndays
integer :: itime

logical :: first_obs

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, prev_time
type(time_type)         :: base_time, delta_time

real(digits12), allocatable :: time(:)
real(r8), allocatable :: lat(:), lon(:)
real(r8), allocatable :: sst(:,:)
real(r8) :: missing_value

!------------------------------------------------------------------------
!  Declare namelist parameters

! sst_error_std ... instrument and representativeness error (std)

real(r8)           :: sst_error_std    = 0.3_r8
character(len=256) :: input_file       = '1234567.nc'
character(len=256) :: output_file_base = 'obs_seq.sst'
logical            :: debug            = .false.
integer            :: subsample_intv   = 1

namelist /oi_sst_to_obs_nml/ input_file, output_file_base, &
                          debug, subsample_intv, sst_error_std

!------------------------------------------------------------------------
! start of executable code
!------------------------------------------------------------------------

! Read and Record the necessary parameters from input.nml
call initialize_utilities()
call find_namelist_in_file('input.nml', 'oi_sst_to_obs_nml', iunit)
read(iunit, nml = oi_sst_to_obs_nml, iostat = io)

if (do_nml_file()) write(nmlfileunit, nml=oi_sst_to_obs_nml)
if (do_nml_term()) write(     *     , nml=oi_sst_to_obs_nml)

call set_calendar_type(GREGORIAN)

ncid  = nc_open_file_readonly(input_file,   routine) 
ndays = nc_get_dimension_size(ncid, 'time', routine)
nlat  = nc_get_dimension_size(ncid, 'lat',  routine)
nlon  = nc_get_dimension_size(ncid, 'lon',  routine)

allocate(time(ndays), lat(nlat), lon(nlon))

call nc_get_variable(ncid, 'time', time, routine) 
call nc_get_variable(ncid, 'lat',  lat,  routine) 
call nc_get_variable(ncid, 'lon',  lon,  routine) 

! ensure longitudes are [0,360] 
where(lon < 0.0_r8) lon = lon + 360.0_r8

base_time = set_base_time(ncid)

num_new_obs = nlon*nlat

allocate(sst(nlon,nlat))

call static_init_obs_sequence()

io = nf90_inq_varid(ncid, 'sst', varid) 
call nc_check(io, routine, context='getting sst variable ID',ncid=ncid)
call nc_get_attribute_from_variable(ncid,'sst','missing_value',missing_value,routine)

TIMELOOP: do itime = 1,ndays

   ! time is stored in the file 2 ways: as real(double) seconds since 1981/1/1,
   ! and as 4 and 2 digit strings for year/mon/day/hr/min
   ! both of these are variables, not attributes

   ! convert to integer days and seconds, and add on to reference time.
   iday = time(itime)
   isec = (time(itime) - iday) * 86400
   delta_time = set_time(isec, iday)
   obs_time = base_time + delta_time
   call get_time(obs_time,  osec, oday)

   ! call print_time(obs_time, str='obs time is ')
   ! call print_date(obs_time, str='obs date is ')

   call get_date(obs_time, year, month, day, hour, minutes, seconds)

   seconds = seconds + (hour*60 + minutes)*60

   write(string1,'(i4,''-'',i2.2,''-'',i2.2,''-'',i5.5)') year,month,day,seconds
   write(output_file,'(A)') trim(output_file_base)//'.'//trim(string1)
   write(*,*)'output file is ',trim(output_file)

   io = nf90_get_var(ncid, varid, sst, start=(/1,1,itime/), count=(/nlon,nlat,1/))
   call nc_check(io, routine, context='get_var sst', ncid=ncid)

   first_obs = .true.
   nmissing = 0
   call init_obs(     obs, num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)
   call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
   call set_copy_meta_data(obs_seq, 1, 'SST observation')
   call   set_qc_meta_data(obs_seq, 1, 'SST QC')

   obslooplat: do j = 1, nlat, subsample_intv
   obslooplon: do i = 1, nlon, subsample_intv
  
     if (sst(i,j) == missing_value) then
        nmissing = nmissing + 1
        cycle obslooplon
     endif
 
     call create_3d_obs(lat(j), lon(i), 0.0_r8, VERTISSURFACE, sst(i,j), &
                         OI_SEA_SURFACE_TEMPERATURE, sst_error_std, oday, osec, qc, obs)
     call add_obs_to_seq(obs_seq, obs, obs_time, prev_obs, prev_time, first_obs)

   enddo obslooplon
   enddo obslooplat

   ! if we added any obs to the sequence, write it out to a file now.
   if ( get_num_obs(obs_seq) > 0 ) then
      if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
      if (debug) print *, '                  skipping = ', nmissing
      call write_obs_seq(obs_seq, output_file)
   else
      write(string1,*)'no observations for output file'
      write(string2,*)'"'//trim(output_file)//'"'
      call error_handler(E_MSG, routine, string1, text2=string2)
   endif

   call destroy_obs_sequence(obs_seq)

end do TIMELOOP

call nc_close_file(ncid,routine)
 
call error_handler(E_MSG, routine, 'Finished successfully.')
call finalize_utilities()

contains


!> 	time:units = "days since 1800-01-01 00:00:00" ;

function set_base_time(ncid)

integer, intent(in) :: ncid
type(time_type)     :: set_base_time

character(len=256) :: timeunits
integer :: io
integer :: year, month, day, hour, minute, second

call nc_get_attribute_from_variable(ncid,'time','units',timeunits,'set_base_time')

read(timeunits,100,iostat=io) year, month, day, hour, minute, second

if (io /= 0) then 
   write(string1,*)'unable to read time base'
   call error_handler(E_ERR, 'set_base_time', string1, &
              source, revision, revdate, text2=timeunits)
endif
            
100 format(11x,i4,5(1x,i2))

set_base_time = set_date(year, month, day, hour, minute, second)

if (debug) then
   write(*,*)'time units is ',trim(timeunits)
   call print_time(set_base_time, str='obs time is ')
   call print_date(set_base_time, str='obs date is ')
endif

end function set_base_time

end program oi_sst_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

!-----------------------------------------------------------------------
!> is2hi_to_obs_netcdf - input is a sea ice thickness file (netcdf format). 
!>                     This program reads the netcdf file and creates an 
!>                     observation sequence file of the sea ice thickness.
!> 
!> Converter credits: Molly M. Wieringa - University of Washington.
!> 
!> Data orginate from https://nsidc.org/data/is2sitdat4/versions/1
!> ICESat-2 L4 Along-Track Sea Ice Thickness, Version 1 (IS2SITDAT4)
!> "This data set reports daily, along-track winter sea ice thickness 
!>  across the Arctic Ocean. Sea ice thickness is estimated using 
!>  ATLAS/ICESat-2 L3A Sea Ice Freeboard (ATL10), Version 5 data and NASA 
!>  Eulerian Snow On Sea Ice Model (NESOSIM) snow loading."
!> 
!> Please make sure you cite the data in accordance to the agreement: 
!>
!> Petty, A. A., N. T. Kurtz, R. Kwok, T. Markus, and T. A. Neumann. 2021. 
!> ICESat-2 L4 Along-Track Sea Ice Thickness, Version 1. [Indicate subset used]. 
!> Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed 
!> Active Archive Center. https://doi.org/10.5067/JTI5YG3S6VAJ. [Date
!> Accessed].
!-----------------------------------------------------------------------

program is2hi_to_obs_netcdf

use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities, finalize_utilities,      &
                              find_namelist_in_file,  &
                              check_namelist_read, nmlfileunit, do_nml_file, &
                              do_nml_term
use  time_manager_mod, only : time_type, set_calendar_type, &
                              set_date, set_time, get_time, GREGORIAN, &
                              operator(>=), operator(-), operator(+)
use      location_mod, only : VERTISSURFACE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,     &
                              static_init_obs_sequence, init_obs,            &
                              write_obs_seq, init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data
use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq, getdimlen
use      obs_kind_mod, only : SAT_SEAICE_AGREG_THICKNESS

use netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, nc_check

use netcdf

implicit none

character(len=*), parameter :: routine = 'is2hi_to_obs_netcdf'

integer :: n, i, j, oday, osec, rcio, iunit, otype, io
integer :: num_copies, num_qc, max_obs, iacc, ialo, ncid, varid
integer :: axdim, aydim
integer :: along_base, across_base
integer, allocatable :: qc_array(:,:)
real(r8), allocatable :: tmask(:,:) ! float tmask:comment = "0 = land, 1 = ocean" ;
character(len=128) :: varname

logical  :: file_exist, first_obs

real(r8) :: temp, qc, wdir, wspeed, werr, thiserr
real(r8) :: uwnd, uerr, vwnd, verr
real(r8) :: dlon, dlat, thislon, thislat
real(r8), allocatable :: lat(:,:), lon(:,:)
real(r8), allocatable :: seaice_thickness(:,:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

! namelist with default values

integer  :: year  = 2018
integer  :: doy   = 287     ! October 14th 
real(r8) :: terr  = 0.10_r8 !MW COMMENT: more presecriptive error possible?
character(len=256) :: seaice_input_file = 'is2data.input'
character(len=256) :: obs_out_file      = 'obs_seq.out'
character(len=256) :: maskfile          = 'cice_hist.nc'
logical  :: debug = .false.  ! set to .true. to print info

namelist /is2hi_to_obs_nc_nml/  year, doy, terr, &
                               seaice_input_file, obs_out_file, &
                               maskfile, debug

! ------------------------
! start of executable code
! ------------------------

call initialize_utilities(routine)

call find_namelist_in_file('input.nml', 'is2hi_to_obs_netcdf_nml', iunit)
read(iunit, nml = is2hi_to_obs_nc_nml, iostat = io)
call check_namelist_read(iunit, io, 'is2hi_to_obs_netcdf_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=is2hi_to_obs_nc_nml)
if (do_nml_term()) write(     *     , nml=is2hi_to_obs_nc_nml)

ncid = nc_open_file_readonly(seaice_input_file, routine)

! get dims along and across the swath path
call getdimlen(ncid, 'ni', axdim)
call getdimlen(ncid, 'nj', aydim)

! remember that when you ncdump the netcdf file, the dimensions are
! listed in C order.  when you allocate them for fortran, reverse the order.
allocate(seaice_thickness(axdim, aydim))
allocate(             lat(axdim, aydim))
allocate(             lon(axdim, aydim))
allocate(           tmask(axdim, aydim))
allocate(        qc_array(axdim, aydim))

varname = 'hi'
io = nf90_inq_varid(ncid, varname, varid)
call nc_check(io, routine, 'nf90_inq_varid "'//trim(varname)//'"')
io = nf90_get_var(ncid, varid, seaice_thickness)
call nc_check(io, routine, 'nf90_get_var "'//trim(varname)//'"')

!! obtain lat and lon
varname = 'lat'
io = nf90_inq_varid(ncid, varname, varid)
call nc_check(io, routine, 'nf90_inq_varid "'//trim(varname)//'"')
io = nf90_get_var(ncid, varid, lat)
call nc_check(io, routine, 'nf90_get_var "'//trim(varname)//'"')

varname = 'lon'
io = nf90_inq_varid(ncid, varname, varid)
call nc_check(io, routine, 'nf90_inq_varid "'//trim(varname)//'"')
io = nf90_get_var(ncid, varid, lon)
call nc_check(io, routine, 'nf90_get_var "'//trim(varname)//'"')

call nc_close_file(ncid, routine, 'data file')

! read in ocean/land mask from a different file.
ncid = nc_open_file_readonly(maskfile, routine)

varname = 'tmask'
io = nf90_inq_varid(ncid, varname, varid)
call nc_check(io, routine, 'nf90_inq_varid "'//trim(varname)//'"')
io = nf90_get_var(ncid, varid, tmask)
call nc_check(io, routine, 'nf90_get_var "'//trim(varname)//'"')

call nc_close_file(ncid, routine, 'mask file')

! convert -180/180 to 0/360
where (lon < 0.0_r8) lon = lon + 360.0_r8

! time setup
call set_calendar_type(GREGORIAN)

!! all obs in a single file are the same time.
comp_day0 = set_date(year, 1, 1, 0, 0, 0)
time_obs = comp_day0 + set_time(0, doy)

! extract time of observation into gregorian day, sec.
call get_time(time_obs, osec, oday)

! There's no actual "vertical" layers. But the 3rd dimension in seaice is category.

! -------------------------

! each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = axdim*aydim
num_copies = 1
num_qc     = 1

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(obs_seq, 1, 'Data QC')

! if you want to append to existing files (e.g. you have a lot of
! small text files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files
! once they are in DART obs_seq format.

!  ! existing file found, append to it
!  inquire(file=obs_out_file, exist=file_exist)
!  if ( file_exist ) then
!     call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)
!  endif

! we have to pick an error range.  since this is a seaice cover fraction
! observation, the valid values should go from 0 to 1.0, so pick 0.1 for now.
qc = 0.0_r8     ! we will reject anything with a bad qc

qc_array = 0    ! making synthetic observations so assume every observation is good

alongloop:  do j = 1, aydim
   acrossloop: do i = 1, axdim

      if (debug) print *, 'start of main loop, ', iacc, ialo

      !! check the lat/lon values to see if they are ok
      if ( lat(i,j) >  90.0_r8 .or. lat(i,j) <   40.0_r8 ) cycle acrossloop
      if ( lon(i,j) <   0.0_r8 .or. lon(i,j) >  360.0_r8 ) cycle acrossloop
      
      ! the actual data values are denser, so inner loop here
      if (qc_array(i,j) /= 0)      cycle acrossloop  !reserve for future quality control
      ! If the mask or data values are outside acceptable bounds, skip them.
      if (   tmask(i,j) <  0.5_r8) cycle acrossloop  !do not convert if it's a land grid
      if (seaice_thickness(i,j) < 0.00_r8) cycle acrossloop
      if (seaice_thickness(i,j) > 25.0_r8) cycle acrossloop

      !>@todo possibly use a higher QC value for suspicious observations
      ! One strategy would be to assign suspicious observations a higher
      ! QC value - this would allow the "input_qc_threshold" namelist to control
      ! whether or not the observation would be assimilated as opposed to having
      ! to modify source code and create multiple versions of the obs_seq file.

      thislat = lat(i,j)
      thislon = lon(i,j)
      ! MW COMMENT: ideally, error for each measurement can be read in from the observation file
      thiserr = terr       ! for now, prescribing a fixed error (instrument) of 10cm
      if (seaice_thickness(i,j) == 0.0_r8) thiserr = 0.05_r8

      ! make an obs derived type, and then add it to the sequence
      call create_3d_obs(thislat, thislon, 0.0_r8, VERTISSURFACE, seaice_thickness(i,j), &
                         SAT_SEAICE_AGREG_THICKNESS, thiserr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      if (debug) print *, 'added seaice obs to output seq'

   end do acrossloop
end do alongloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

deallocate(seaice_thickness, lon, lat, qc_array, tmask)

! end of main program
call finalize_utilities()

end program is2hi_to_obs_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

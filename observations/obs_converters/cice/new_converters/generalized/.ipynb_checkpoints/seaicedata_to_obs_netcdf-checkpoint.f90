! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

!-----------------------------------------------------------------------
!> seaicedata_to_obs_netcdf - input is a sea ice data file in netcdf 
!>                            format that has been generated using {INSERT
!>                            NEW TOOL HERE}. This program reads the netcdf 
!>                            file and creates an observation sequence file
!>                            of the sea ice data.
!> 
!> Converter credits: Molly M. Wieringa - University of Washington.
!> 
!> This data converter is designed to be general, with the data type, 
!> values, and errors being read in from the data file or namelist. The  
!> data file should be generated using the associated python program, which 
!> standardizes the data and its naming conventions and metadata. 
!-----------------------------------------------------------------------

program seaicedata_to_obs_netcdf

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
use      obs_kind_mod, only : get_name_for_type_of_obs

use netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, nc_check

use netcdf

implicit none

character(len=*), parameter :: routine = 'seaicedata_to_obs_netcdf'

integer :: n, i, j, oday, osec, rcio, iunit, otype, io
integer :: num_copies, num_qc, max_obs, iacc, ialo, ncid, varid, data_type_int
integer :: axdim, aydim, catdim, len_time
integer :: along_base, across_base
real(r8), allocatable :: tmask(:,:) ! float tmask:comment = "0 = land, 1 = ocean" ;
character(len=128) :: varname

logical  :: file_exist, first_obs

real(r8) :: temp, qc, wdir, wspeed, werr, thiserr
real(r8) :: uwnd, uerr, vwnd, verr
real(r8) :: dlon, dlat, thislon, thislat
real(r8), allocatable :: lat(:,:,:), lon(:,:,:)
real(r8), allocatable :: seaiceerror(:,:,:,:)
real(r8), allocatable :: seaicedata(:,:,:,:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

! namelist with default values
integer  :: year      = 2000 ! beginning of desired reconstruction
integer  :: doy       = 1    ! October 14th 
integer  :: data_type = 12   ! integer index of desired data type, SIC by default
character(len=256) :: seaice_input_file = 'seaicedata.input'
character(len=256) :: obs_out_file      = 'obs_seq.out'
character(len=256) :: maskfile          = 'cice_hist.nc'
logical  :: itd_ob = .false. ! set to .true. if observation has ITD dimension
logical  :: debug = .false.  ! set to .true. to print info

namelist /seaicedata_to_obs_netcdf_nml/  year, doy, data_type, seaice_input_file, &
                                         obs_out_file, maskfile, itd_ob, &
                                         debug

! ------------------------
! start of executable code
! ------------------------

call initialize_utilities(routine)

call find_namelist_in_file('input.nml', 'seaicedata_to_obs_netcdf_nml', iunit)
read(iunit, nml = seaicedata_to_obs_netcdf_nml, iostat = io)
call check_namelist_read(iunit, io, 'seaicedata_to_obs_netcdf_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=seaicedata_to_obs_netcdf_nml)
if (do_nml_term()) write(     *     , nml=seaicedata_to_obs_netcdf_nml)

ncid = nc_open_file_readonly(seaice_input_file, routine)

! get dims along and across the swath path
call getdimlen(ncid, 'time', len_time)
call getdimlen(ncid,  'lat',    axdim)
call getdimlen(ncid,  'lon',    aydim)
if (itd_ob) then 
   call getdimlen(ncid, 'nc', catdim)
else 
   catdim = 1
endif

! remember that when you ncdump the netcdf file, the dimensions are
! listed in C order.  when you allocate them for fortran, reverse the order.
allocate( seaice_data(len_time, axdim, aydim, catdim))
allocate(seaice_error(len_time, axdim, aydim, catdim))
allocate(         lat(len_time, axdim, aydim))
allocate(         lon(len_time, axdim, aydim))
allocate(       tmask(axdim, aydim))


! get string name type of observation
OBS_TYPE = get_name_for_type_of_obs(data_type)

! get data
varname = 'data'
io = nf90_inq_varid(ncid, varname, varid)
call nc_check(io, routine, 'nf90_inq_varid "'//trim(varname)//'"')
io = nf90_get_var(ncid, varid, seaice_data)
call nc_check(io, routine, 'nf90_get_var "'//trim(varname)//'"')

! get errors
varname = 'error'
io = nf90_inq_varid(ncid, varname, varid)
call nc_check(io, routine, 'nf90_inq_varid "'//trim(varname)//'"')
io = nf90_get_var(ncid, varid, seaice_error)
call nc_check(io, routine, 'nf90_get_var "'//trim(varname)//'"')

! get lat and lon
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

! close file
call nc_close_file(ncid, routine, 'data file')

! read in ocean/land mask from a different file.
ncid = nc_open_file_readonly(maskfile, routine)

varname = 'tmask'
io = nf90_inq_varid(ncid, varname, varid)
call nc_check(io, routine, 'nf90_inq_varid "'//trim(varname)//'"')
io = nf90_get_var(ncid, varid, tmask)
call nc_check(io, routine, 'nf90_get_var "'//trim(varname)//'"')

! close mask file
call nc_close_file(ncid, routine, 'mask file')

! convert -180/180 to 0/360
where (lon < 0.0_r8) lon = lon + 360.0_r8

! time setup
call set_calendar_type(GREGORIAN)

! -------------------------

! each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = axdim*aydim*catdim
num_copies = 1
num_qc     = 1

do t = 1, len_time

   ! all obs in a single file are the same time.
   comp_day0 = set_date(year, 1, 1, 0, 0, 0)
   time_obs = comp_day0 + set_time(0, doy)

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)
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

   ! we have to pick an error range.  since this is a seaice cover fraction
   ! observation, the valid values should go from 0 to 1.0, so pick 0.1 for now.
   qc = 0.0_r8     ! we will reject anything with a bad qc

   ! move through the observations and create a DART 3d observation if they pass
   ! quality control checks 
   alongloop:  do j = 1, aydim
      acrossloop: do i = 1, axdim
         do k = 1, catdim

            if (debug) print *, 'start of main loop, ', iacc, ialo

            !! check the lat/lon values to see if they are ok
            if ( lat(t,i,j) >  90.0_r8 .or. lat(t,i,j) <   40.0_r8 ) cycle acrossloop
            if ( lon(t,i,j) <   0.0_r8 .or. lon(t,i,j) >  360.0_r8 ) cycle acrossloop
      
            ! If the mask or data values are outside acceptable bounds, skip them.
            if (   tmask(i,j) <  0.5_r8) cycle acrossloop  !do not convert if it's a land grid
            if ( seaice_data(t,i,j,k) < 0.00_r8 .or.  seaice_data(t,i,j,k) > 25.0_r8) cycle acrossloop
            if (seaice_error(t,i,j,k) < 0.00_r8 ) cycle acrossloop

            ! assign latitude, longitude, category, and error of okay values
            thislat = lat(t,i,j)
            thislon = lon(t,i,j)
            thiscat = k - 1 
            thiserr = seaice_error(t,i,j,k)     

            ! make an obs derived type, and then add it to the sequence
            call create_3d_obs(thislat, thislon, thiscat, VERTISSURFACE, seaice_data(t,i,j,k), &
                               data_type, thiserr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

            if (debug) print *, 'added ', OBS_TYPE, ' obs to output seq'
         end do
      end do acrossloop
   end do alongloop

   ! if we added any obs to the sequence, write it out to a file now.
   if ( get_num_obs(obs_seq) > 0 ) then
      if (debug) print *, 'writing obs_seq for ', time_obs, ', obs_count = ', get_num_obs(obs_seq)
      call write_obs_seq(obs_seq, obs_out_file)
   endif

   ! adjust time to move through next timestep
   doy = doy + 1
   if (doy > 365) then
      doy = 1
      year = year + 1
   end if

enddo 

! release allocated arrays 
deallocate(seaice_data, seaice_error, lon, lat, tmask)

! end of main program
call finalize_utilities()

end program seaicedata_to_obs_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$


program seaice_to_obs_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   seaice_to_obs_netcdf - input is a seaice-coverage file that has been
!      converted from HDF to netCDF with an automated tool.  this
!      program then takes the unsigned byte/integer(1) data and 
!      converts it into a seaice coverage obs_seq file.
!
!     created  5 jul 2012   nancy collins NCAR/IMAGe
!     based on a previous text-based version done in mar 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD
use     utilities_mod, only : initialize_utilities, finalize_utilities,      &
                              open_file, close_file, find_namelist_in_file,  &
                              check_namelist_read, nmlfileunit, do_nml_file, &
                              do_nml_term, nc_check
use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date
use      location_mod, only : VERTISLEVEL
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,     &
                              static_init_obs_sequence, init_obs,            &
                              write_obs_seq, init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq, getdimlen
use      obs_kind_mod, only : SYN_SEAICE_CONCENTR

use netcdf

implicit none

!character(len=64), parameter :: obs_out_file    = 'obs_seq.out'

integer :: n, i, j, oday, osec, rcio, iunit, otype, io
integer :: num_copies, num_qc, max_obs, iacc, ialo, ncid, varid
integer :: axdim, aydim
integer :: along_base, across_base
integer, allocatable :: qc_array(:,:)
integer, allocatable :: tmask(:,:)           
character(len=128) :: varname

logical  :: file_exist, first_obs

real(r8) :: temp, qc, wdir, wspeed, werr
real(r8) :: uwnd, uerr, vwnd, verr
real(r8) :: dlon, dlat, thislon, thislat
real(r8), allocatable :: lat(:,:), lon(:,:)
real(r8), allocatable :: seaice_concentr(:,:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

integer  :: year  = 2000
integer  :: doy   = 1
real(r8) :: terr = 0.1_r8   ! FIXME - wild guess
real(r8) :: cat  = 1_r8
character(len=128) :: seaice_input_file = 'seaicedata.input'
character(len=128) :: obs_out_file    = 'obs_seq.out'
character(len=128) :: maskfile        = 'cice_hist.nc'
logical  :: debug = .false.  ! set to .true. to print info


namelist /seaice_to_obs_nc_nml/  year, doy, cat, terr, &
                               seaice_input_file, obs_out_file, &
                               maskfile, debug

! ------------------------
! start of executable code
! ------------------------

call initialize_utilities('seaice_to_obs_netcdf')

call find_namelist_in_file('input.nml', 'seaice_to_obs_nc_nml', iunit)
read(iunit, nml = seaice_to_obs_nc_nml, iostat = io)
call check_namelist_read(iunit, io, 'seaice_to_obs_nc_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=seaice_to_obs_nc_nml)
if (do_nml_term()) write(     *     , nml=seaice_to_obs_nc_nml)

! open netcdf file here.
call nc_check( nf90_open(seaice_input_file, nf90_nowrite, ncid), &
               'seaice_to_obs_netcdf', 'opening file '//trim(seaice_input_file))

! get dims along the swath path, and across the swath path.  the rest of
! the data arrays use these for their dimensions
call getdimlen(ncid, 'ni', axdim)
call getdimlen(ncid, 'nj', aydim)

! remember that when you ncdump the netcdf file, the dimensions are
! listed in C order.  when you allocate them for fortran, reverse the order.
allocate(seaice_concentr(axdim, aydim))
allocate(lon(axdim,aydim), lat(axdim,aydim))
allocate(qc_array(axdim,aydim))
allocate(tmask(axdim,aydim))

varname = 'aice_obs'
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'seaice_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, seaice_concentr), &
               'seaice_to_obs_netcdf', 'getting var '// trim(varname))

!! obtain lat and lon 
varname = 'lat'
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'seaice_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, lat), &
               'seaice_to_obs_netcdf', 'getting var '// trim(varname))

varname = 'lon'
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'seaice_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, lon), &
               'seaice_to_obs_netcdf', 'getting var '// trim(varname))

! read in ocean/land mask
varname = 'tmask'
call nc_check( nf90_open(maskfile,nf90_nowrite,ncid), &
               'seaice_to_obs_netcdf','opening file'//trim(maskfile))
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'seaice_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, tmask), &
               'seaice_to_obs_netcdf', 'getting var '// trim(varname))

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
      if ( lat(i,j) >  90.0_r8 .or. lat(i,j) <  40.0_r8 ) cycle alongloop
      if ( lon(i,j) <   0.0_r8 .or. lon(i,j) >  360.0_r8 ) cycle acrossloop
    
      ! the actual data values are denser, so inner loop here
            
            if (qc_array(i,j) /= 0) cycle alongloop  !reserve for future quality control
            if (tmask(i,j)/= 1) cycle alongloop  !do not convert if it's a land grid
            if (seaice_concentr(i,j).lt.0.01_r8) cycle acrossloop   !FIXME temporary do not assimilate 
                                                                    !when observed sea ice is 0 coverage
            ! compute the lat/lon for this obs  FIXME: this isn't right

            thislat = lat(i,j)

            thislon = lon(i,j)

            ! make an obs derived type, and then add it to the sequence
            call create_3d_obs(thislat, thislon, cat, VERTISLEVEL, seaice_concentr(i,j), &
                               SYN_SEAICE_CONCENTR, terr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
         
            if (debug) print *, 'added seaice obs to output seq'

   end do acrossloop
end do alongloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

end program seaice_to_obs_netcdf

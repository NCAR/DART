! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program snow_to_obs_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   snow_to_obs_netcdf - input is a snow-coverage file that has been
!      converted from HDF to netCDF with an automated tool.  this
!      program then takes the unsigned byte/integer(1) data and 
!      converts it into a snow coverage obs_seq file.
!
!     created  5 jul 2012   nancy collins NCAR/IMAGe
!     based on a previous text-based version done in mar 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD
use     utilities_mod, only : initialize_utilities, finalize_utilities,      &
                              open_file, close_file, find_namelist_in_file,  &
                              check_namelist_read, nmlfileunit, do_nml_file, &
                              do_nml_term
use  netcdf_utilities_mod, only : nc_check
use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date
use      location_mod, only : VERTISSURFACE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,     &
                              static_init_obs_sequence, init_obs,            &
                              write_obs_seq, init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq, getdimlen
use      obs_kind_mod, only : MODIS_SNOWCOVER_FRAC

use netcdf

implicit none

character(len=64), parameter :: obs_out_file    = 'obs_seq.out'

integer :: n, i, j, oday, osec, rcio, iunit, otype, io
integer :: num_copies, num_qc, max_obs, iacc, ialo, ncid, varid
integer :: along, across, coarse_along, coarse_across, per_along, per_across
integer :: along_base, across_base
integer, allocatable :: qc_array(:,:), snowcover_type(:,:)
           
integer(1), allocatable :: sn_byte(:,:)
character(len=128) :: varname

logical  :: file_exist, first_obs

real(r8) :: temp, terr, qc, wdir, wspeed, werr
real(r8) :: vert, uwnd, uerr, vwnd, verr
real(r8) :: dlon, dlat, thislon, thislat
real(r8), allocatable :: lat(:,:), lon(:,:)
real(r8), allocatable :: fsnowcover(:,:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

integer  :: year  = 2000
integer  :: doy   = 1
character(len=128) :: snow_input_file = 'snowdata.input'
logical  :: debug = .false.  ! set to .true. to print info


namelist /snow_to_obs_nc_nml/  year, doy, &
                               snow_input_file, debug

! ------------------------
! start of executable code
! ------------------------

call initialize_utilities('snow_to_obs_netcdf')

call find_namelist_in_file('input.nml', 'snow_to_obs_nc_nml', iunit)
read(iunit, nml = snow_to_obs_nc_nml, iostat = io)
call check_namelist_read(iunit, io, 'snow_to_obs_nc_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=snow_to_obs_nc_nml)
if (do_nml_term()) write(     *     , nml=snow_to_obs_nc_nml)

! open netcdf file here.
call nc_check( nf90_open(snow_input_file, nf90_nowrite, ncid), &
               'snow_to_obs_netcdf', 'opening file '//trim(snow_input_file))

! get dims along the swath path, and across the swath path.  the rest of
! the data arrays use these for their dimensions
call getdimlen(ncid, 'Along_swath_lines_500m_MOD_Swath_Snow', along)
call getdimlen(ncid, 'Cross_swath_pixels_500m_MOD_Swath_Snow', across)
call getdimlen(ncid, 'Coarse_swath_lines_5km_MOD_Swath_Snow', coarse_along)
call getdimlen(ncid, 'Coarse_swath_pixels_5km_MOD_Swath_Snow', coarse_across)

! remember that when you ncdump the netcdf file, the dimensions are
! listed in C order.  when you allocate them for fortran, reverse the order.
allocate(sn_byte(across, along))
allocate(fsnowcover(across, along))
allocate(snowcover_type(across, along))
allocate(qc_array(across, along))
allocate(lon(coarse_across, coarse_along), lat(coarse_across, coarse_along))

! apparently the lat/lon info is every Nth data value in each dim
! FIXME: it seems like there are 2 less across values than there
! should be to have 10 obs / across swath??
per_along = along / coarse_along
per_across = (across+2) / coarse_across

! read the data for the requested arrays.
varname = 'Fractional_Snow_Cover'
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'snow_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, sn_byte), &
               'snow_to_obs_netcdf', 'getting var '// trim(varname))

! snowcover fraction is stored as unsigned bytes in the netcdf file.
! read into a 1 byte int array, convert to unsigned and real
! values from 0 to 100 are valid percentages, anything from 101 to 256
! is a flag (see the netcdf header for details) and will be ignored.

fsnowcover = sn_byte
where (fsnowcover < 0.0_r8) fsnowcover = fsnowcover + 256.0_r8

varname = 'Snow_Cover'
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'snow_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, sn_byte), &
               'snow_to_obs_netcdf', 'getting var '// trim(varname))

varname = 'Snow_Cover_Pixel_QA'
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'snow_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, sn_byte), &
               'snow_to_obs_netcdf', 'getting var '// trim(varname))

qc_array = sn_byte
where (qc_array < 0) qc_array = qc_array + 256

varname = 'Longitude'
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'snow_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, lon), &
               'snow_to_obs_netcdf', 'getting var '// trim(varname))

! convert -180/180 to 0/360
where (lon < 0.0_r8) lon = lon + 180.0_r8

varname = 'Latitude'
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'snow_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, lat), &
               'snow_to_obs_netcdf', 'getting var '// trim(varname))

! not sure if we need this one but read it in just in case
! this seems to be some kind of cell type (night, no snow, lake, ocean,
! cloud, lake ice snow, detector saturated, missing data, fill)
varname = 'Snow_Cover'
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'snow_to_obs_netcdf', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, sn_byte), &
               'snow_to_obs_netcdf', 'getting var '// trim(varname))

! snowcover is stored as unsigned bytes in the netcdf file.
! read into a character array and convert to integer.

snowcover_type = sn_byte
where (snowcover_type < 0) snowcover_type = snowcover_type + 256


! time setup
call set_calendar_type(GREGORIAN)

!! all obs in a single file are the same time.
comp_day0 = set_date(year, 1, 1, 0, 0, 0)
time_obs = comp_day0 + set_time(0, doy)

! extract time of observation into gregorian day, sec.
call get_time(time_obs, osec, oday)

! surface obs.  normally we set vert to the actual surface elevation, 
! we do not have it here, so set to 0 for now.
vert = 0.0_r8
   
! -------------------------

! each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = across * along
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

! we have to pick an error range.  since this is a snow cover fraction
! observation, the valid values should go from 0 to 1.0, so pick 0.1 for now.
terr = 0.1_r8   ! FIXME - wild guess
qc = 0.0_r8     ! we will reject anything with a bad qc


alongloop:  do ialo = 1, coarse_along

   acrossloop: do iacc = 1, coarse_across

if (debug) print *, 'start of main loop, ', iacc, ialo

      !! check the lat/lon values to see if they are ok
      if ( lat(iacc, ialo) >  90.0_r8 .or. lat(iacc, ialo) <  -90.0_r8 ) cycle alongloop
      if ( lon(iacc, ialo) <   0.0_r8 .or. lon(iacc, ialo) >  360.0_r8 ) cycle alongloop
    
      ! the actual data values are denser, so inner loop here
      along_base = (ialo-1) * per_along + 1
      inner_along: do j = along_base, along_base + per_along
 
         across_base = (iacc-1) * per_across + 1
         inner_across: do i = across_base, across_base + per_across

            if (i > across) cycle inner_across
            if (j > along) cycle inner_across
            if (qc_array(i,j) /= 0) cycle inner_across
            if (fsnowcover(i, j) > 101) cycle inner_across

            ! compute the lat/lon for this obs  FIXME: this isn't right

            thislat = lat(iacc, ialo)
            thislon = lon(iacc, ialo)

            ! make an obs derived type, and then add it to the sequence
            call create_3d_obs(thislat, thislon, vert, VERTISSURFACE, fsnowcover(i,j)/100.0_r8, &
                               MODIS_SNOWCOVER_FRAC, terr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
         
            if (debug) print *, 'added snow obs to output seq'

         end do inner_across
      end do inner_along
   end do acrossloop
end do alongloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

end program snow_to_obs_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

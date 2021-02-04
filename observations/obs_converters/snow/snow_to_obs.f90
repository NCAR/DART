! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program snow_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   snow_to_obs - a program that only needs minor customization to read
!      in a text-based dataset - either white-space separated values or
!      fixed-width column data.
!
!     created 29 Mar 2010   nancy collins NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD
use     utilities_mod, only : initialize_utilities, finalize_utilities,      &
                              open_file, close_file, find_namelist_in_file,  &
                              check_namelist_read, nmlfileunit, do_nml_file, &
                              do_nml_term
use  time_manager_mod, only : time_type, set_calendar_type, set_date, set_time, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date
use      location_mod, only : VERTISSURFACE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,     &
                              static_init_obs_sequence, init_obs,            &
                              write_obs_seq, init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq
use      obs_kind_mod, only : MODIS_SNOWCOVER_FRAC

implicit none

character(len=64), parameter :: obs_out_file    = 'obs_seq.out'

integer :: oday, osec, rcio, iunit, io
integer :: num_copies, num_qc, max_obs, ix, iy
           
logical  :: file_exist, first_obs

real(r8) :: terr, qc
real(r8) :: lat, lon, vert
real(r8) :: dlon, dlat
real(r8), allocatable :: coverage(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

integer  :: longrid = 360
integer  :: latgrid = 90
integer  :: year  = 2000
integer  :: doy   = 1
character(len=128) :: snow_input_file = 'snowdata.input'
real(r8) :: missing_value = -20.0_r8
logical  :: debug = .false.  ! set to .true. to print info


namelist /snow_to_obs_nml/  longrid, latgrid, year, doy, &
                            snow_input_file, missing_value, debug

! ------------------------
! start of executable code
! ------------------------

call initialize_utilities('snow_to_obs')

call find_namelist_in_file('input.nml', 'snow_to_obs_nml', iunit)
read(iunit, nml = snow_to_obs_nml, iostat = io)
call check_namelist_read(iunit, io, 'snow_to_obs_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=snow_to_obs_nml)
if (do_nml_term()) write(     *     , nml=snow_to_obs_nml)


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
   
! make space for values
allocate(coverage(longrid))

! -------------------------

! open input snow file

iunit = open_file(snow_input_file, 'formatted', 'read')
if (debug) print *, 'opened input file ' // trim(snow_input_file)


! each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = 100000
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

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8
terr = 1.0_r8   ! FIXME - wild guess

! just northern hemisphere
dlon = 360.0_r8 / longrid
dlat =  90.0_r8 / latgrid


lat = 0.0_r8
latloop: do iy = 1, latgrid

if (debug) print *, 'start of latloop, iy = ', iy
   read(iunit, *, iostat=rcio) coverage

   if (rcio /= 0) then 
      if (debug) print *, 'got bad read code getting line, rcio = ', rcio
      exit latloop
   endif
   

   lon = 0.0_r8

   lonloop:  do ix = 1, longrid


      !! check the lat/lon values to see if they are ok
      !if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle lonloop
      !if ( lon <   0.0_r8 .or. lon >  360.0_r8 ) cycle lonloop
   
      if(coverage(ix) == missing_value) cycle lonloop

if (debug) print *, ix, iy, 'got coverage ', coverage(ix)

      ! make an obs derived type, and then add it to the sequence
      call create_3d_obs(lat, lon, vert, VERTISSURFACE, coverage(ix), &
                         MODIS_SNOWCOVER_FRAC, terr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
      if (debug) print *, 'added snow obs to output seq'

      lon = lon + dlon
   end do lonloop
if (debug) print *, 'end of lonloop'
   
   lat = lat + dlat

end do latloop
if (debug) print *, 'end of latloop'

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

end program snow_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

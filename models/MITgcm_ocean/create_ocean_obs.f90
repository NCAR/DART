! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program create_ocean_obs

! Initial program to read the raw ocean observations and insert them
! into an observation sequence. To make things easy ... we will mandate
! an assimilation interval of 1 day - so all observations will be
! redefined to occur at NOON on the day they were observed.

use         types_mod, only : r8, deg2rad, PI

use  obs_sequence_mod, only : obs_sequence_type, write_obs_seq, &
                              static_init_obs_sequence, destroy_obs_sequence 

use dart_MITocean_mod, only : real_obs_sequence

use     utilities_mod, only : initialize_utilities, register_module, &
                              do_output, logfileunit, &
                              error_handler, finalize_utilities, E_ERR, E_MSG, &
                              find_namelist_in_file, check_namelist_read

use  time_manager_mod, only : time_type, set_date, set_time, print_date, &
                              operator(+), set_calendar_type, GREGORIAN
  

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'create_ocean_obs.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

type(obs_sequence_type) :: seq

integer :: iunit, io
type(time_type) :: time1, timeN

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------
        
integer :: year = 1996, month = 1, day = 1, tot_days = 31
integer :: max_num = 800000
character(len=256) :: fname = 'raw_ocean_obs.txt'
character(len=256) :: output_name = 'raw_ocean_obs_seq.out'
logical :: hfradar = .false.

real(r8) :: lon1 =   0.0_r8,  &   !  lower longitude bound
            lon2 = 360.0_r8,  &   !  upper longitude bound 
            lat1 = -90.0_r8,  &   !  lower latitude bound
            lat2 =  90.0_r8       !  upper latitude bound

namelist /create_ocean_obs_nml/ year, month, day, tot_days, max_num, &
        fname, output_name, lon1, lon2, lat1, lat2, hfradar

! ----------------------------------------------------------------------
! start of executable program code
! ----------------------------------------------------------------------

call initialize_utilities('create_ocean_obs')
call register_module(source,revision,revdate)

! Initialize the obs_sequence module ...

call static_init_obs_sequence()

! Read the namelist entry
call find_namelist_in_file("input.nml", "create_ocean_obs_nml", iunit)
read(iunit, nml = create_ocean_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "create_ocean_obs_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'create_ocean_obs','create_ocean_obs_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=create_ocean_obs_nml)
if (do_output()) write(     *     , nml=create_ocean_obs_nml)

! Real observations are required to use the Gregorian calendar.
call set_calendar_type(GREGORIAN)
time1 = set_date(year, month, day)
timeN = time1 + set_time(0,tot_days)
call print_date(time1,str='First date of interest')
call print_date(timeN,str='Last  date of interest')

! The file is read and parsed into a DART observation sequence linked list
seq = real_obs_sequence(fname, time1, timeN, max_num, &
                         lon1, lon2, lat1, lat2, hfradar)

call write_obs_seq(seq, output_name)

call destroy_obs_sequence(seq) ! release the memory of the seq.

call error_handler(E_MSG,'create_ocean_obs','Finished successfully.',source,revision,revdate)
call finalize_utilities()

end program create_ocean_obs


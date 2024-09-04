! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! This file is meant to read a text file containing bottle data from the
! Bermuda Atlantic Time-Series Study (https://bats.bios.asu.edu/), which
! is converted to climatological estimates with monthly resolution.

program bats_to_clim_obs

use         types_mod, only : r8
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file, &
                              find_namelist_in_file, check_namelist_read, &
                              error_handler, E_ERR, E_MSG, nmlfileunit,   &
                              do_nml_file, do_nml_term
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), &
                              print_date
use      location_mod, only : VERTISHEIGHT, VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, set_copy_meta_data, &
                              set_qc_meta_data, destroy_obs_sequence
use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq
use      obs_kind_mod, only : BATS_OXYGEN, BATS_INORGANIC_CARBON, BATS_ALKALINITY, &
                              BATS_NITRATE, BATS_PHOSPHATE, BATS_SILICATE

implicit none

integer, parameter :: NUM_SCALAR_OBS = 6  ! maximum number of scalar observation variables that will
                                          ! be assimilated at each observation.

! this array defines the order in which observations are read from the file
integer, parameter :: OTYPE_ORDERING(NUM_SCALAR_OBS) &
                      = (/BATS_OXYGEN, BATS_INORGANIC_CARBON, BATS_ALKALINITY, BATS_NITRATE, &
                          BATS_PHOSPHATE, BATS_SILICATE/)

real(r8), parameter :: MIN_OBS_ERROR = 0.1_r8
real(r8), parameter :: BATS_LON      = 360.0_r8 - 64.0_r8
real(r8), parameter :: BATS_LAT      = 31.0_r8

! namelist variables, changeable at runtime
character(len=256) :: text_input_file       = 'bats_climatology.txt'
character(len=256) :: obs_out_dir           = 'obs_seq_files'
integer            :: max_lines             = 3000
real(r8)           :: obs_err_var_inflation = 10.0
logical            :: debug                 = .true.

namelist /bats_to_clim_obs_nml/ text_input_file, max_lines, obs_err_var_inflation, obs_out_dir, debug

! local variables
character (len=294) :: input_line, obs_out_file
character (len=6)   :: month_str

integer  :: oday, osec, rcio, iunit, char_index, comma_index, line_number, otype_index, &
            month, month_old, num_copies, num_qc, max_obs
integer  :: comma_locs(4)
logical  :: first_obs, new_obs_seq

real(r8) :: qc, obs_err
real(r8) :: vert, ovalue

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time


! start of executable code

call initialize_utilities('bats_to_clim_obs')

! Read the namelist entries
call find_namelist_in_file("input.nml", "bats_to_clim_obs_nml", iunit)
read(iunit, nml = bats_to_clim_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "bats_to_clim_obs_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=bats_to_clim_obs_nml)
if (do_nml_term()) write(     *     , nml=bats_to_clim_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)

! open input text file

iunit = open_file(text_input_file, 'formatted', 'read')
if (debug) print *, 'opened input file ' // trim(text_input_file)

max_obs    = NUM_SCALAR_OBS*max_lines
num_copies = 1
num_qc     = 1

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! Set the DART data quality control.
qc = 0.0_r8

! used later to manage splitting of data into obs-sequence files
month_old = 0

line_number = 0 ! counts the number of lines that have been read so far

obsloop: do    ! no end limit - have the loop break when input ends
   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=rcio) input_line
   line_number = line_number + 1

   if (rcio /= 0) then 
      if (debug) print *, 'got bad read code from input file at line ',line_number,', rcio = ', rcio
      exit obsloop
   endif

   ! finding the locations of commas within the string

   comma_index = 0
   char_index  = 0

   do while(comma_index < 4)
      char_index = char_index + 1

      if(input_line(char_index:char_index) ==  ',') then
         comma_index = comma_index + 1
         comma_locs(comma_index) = char_index
      end if
   end do

   ! reading the month value of this line

   read(input_line(1:(comma_locs(1) - 1)), *, iostat=rcio) month
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting month at line ",line_number,", rcio = ",rcio
      exit obsloop
   end if

   ! reading the observation type

   read(input_line((comma_locs(1) + 1):(comma_locs(2) - 1)), *, iostat=rcio) otype_index
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting obs type at line ",line_number,", rcio = ",rcio
      exit obsloop
   end if

   ! reading the observation depth

   read(input_line((comma_locs(2) + 1):(comma_locs(3) - 1)), *, iostat=rcio) vert
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting vert at line ",line_number,", rcio = ",rcio
      exit obsloop
   end if

   ! reading the observtaion value

   read(input_line((comma_locs(3) + 1):(comma_locs(4) - 1)), *, iostat=rcio) ovalue
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting observation value at line ",line_number,", rcio = ",rcio
      exit obsloop
   end if

   ! reading the observtaion uncertainty (standard deviation)

   read(input_line((comma_locs(4) + 1):len(input_line)), *, iostat=rcio) obs_err
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting observation uncertainty at line ",line_number,", rcio = ",rcio
      exit obsloop
   end if

   ! inflating the observation error according to the number of cycles in MARBL-DART

   obs_err = obs_err*sqrt(obs_err_var_inflation)

   ! unit corrections from BATS to MARBL

   if((OTYPE_ORDERING(otype_index) == BATS_ALKALINITY) .or. (OTYPE_ORDERING(otype_index) == BATS_INORGANIC_CARBON)) then
      ovalue  = ovalue*1.026
      obs_err = obs_err*1.026
   end if

   ! For ensemble smoothing, each observation is considered to be a measurement
   ! of the same climatological "year", which is arbitrarily labeled as year 1601
   ! (the first year in DART's internal calendar). Observations are given identical
   ! day/hour/minute timestamps to avoid time-ordering errors in the obs-seq file.

   time_obs = set_date(1601, 1, 1, hours=0, minutes=0)
   call get_time(time_obs, osec, days=oday)

   if(debug) then
      print *, "adding climatological mean for month ",month
      print *, " \__ type:        ",OTYPE_ORDERING(otype_index)
      print *, "     vert:        ",vert
      print *, "     value:       ",ovalue
      print *, "     uncertainty: ",obs_err
      print *, "     oday:        ",oday
      print *, "     osec:        ",osec
   end if

   new_obs_seq = (first_obs .or. (month /= month_old))
   month_old   = month

   ! if necessary, saving the old observation sequence and beginning a new one.
   if(new_obs_seq) then
      ! dumping the observations so far into their own file
      if ( (.not. first_obs) .and. (get_num_obs(obs_seq) > 0) ) then
         if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
         call write_obs_seq(obs_seq, obs_out_file)
         call destroy_obs_sequence(obs_seq)
         first_obs = .true.
      endif

      write(month_str, "(I0.2)") month
      obs_out_file = trim(obs_out_dir)//"/BATS_clim_"//trim(month_str)//".out"

      ! create a new, empty obs_seq file.
      call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

      ! the first one needs to contain the string 'observation' and the
      ! second needs the string 'QC'.
      call set_copy_meta_data(obs_seq, 1, 'observation')
      call set_qc_meta_data(obs_seq, 1, 'Data QC')
   end if

   ! adding the observation
   call create_3d_obs(BATS_LAT, BATS_LON, vert, VERTISHEIGHT, &
                         ovalue, OTYPE_ORDERING(otype_index), obs_err, &
                         oday, osec, qc, obs)

   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
end do obsloop

! putting any remaining observations into an obs sequence file

if ( (.not. first_obs) .and. (get_num_obs(obs_seq) > 0) ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
   call destroy_obs_sequence(obs_seq)
endif

! end of main program
call finalize_utilities()

end program bats_to_clim_obs

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program ascii_table_converter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   ascii_table_converter - a program that only needs minor customization
!      to read in a text-based dataset - either white-space separated 
!      values or fixed-width column data.
!
! time latitude longitude soil_moisture soil_moisture_error retrieval_quality_flag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, digits12

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              register_module, error_handler, E_MSG, E_ERR, &
                              open_file, close_file, do_nml_file, do_nml_term, &
                              check_namelist_read, find_namelist_in_file, &
                              nmlfileunit, logfileunit, file_exist, find_textfile_dims

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, set_time, get_time, print_time, &
                              print_date, operator(-), operator(+), operator(>), &
                              operator(<), operator(==), operator(<=), operator(>=)

use      location_mod, only : VERTISHEIGHT

use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : SOIL_MOISTURE

implicit none

character(len=*), parameter :: source   = 'ascii_table_converter.f90'
character(len=*), parameter :: routine='ascii_table_converter'

!-----------------------------------------------------------------------
! Namelist with default values
!-----------------------------------------------------------------------

character(len=256) :: input_file_list = 'file_list.txt'
character(len=256) :: obs_out_file = 'obs_seq.out'

real(r8) :: maxgoodqc = 3.0_r8
logical  :: verbose   = .true.

namelist /ascii_table_converter_nml/ &
         input_file_list, obs_out_file, maxgoodqc, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

integer, parameter :: MAX_NUM_INPUT_FILES = 500 ! max number of input files to be processed
integer, parameter :: NOBS_PER_FILE = 20000     ! max number of observations per file

character(len=256), dimension(MAX_NUM_INPUT_FILES) :: filename_seq_list
character(len=256) :: filename
character(len=512) :: string1, string2, string3

integer  :: num_input_files = 0  ! actual number of files
integer  :: i, ifile, iline
integer  :: num_obs_used, num_obs_read
integer  :: days, seconds, iocode, iunit
integer  :: num_copies, num_qc, max_obs

logical  :: first_obs

type(obs_sequence_type) :: obs_seq

type(obs_type)  :: obs, prev_obs

type(time_type) :: reference_time, offset, obs_time, prev_time

real(digits12) :: timestamp
real(r8) :: latitude
real(r8) :: longitude
real(r8) :: observation
real(r8) :: soil_moisture_error
real(r8) :: qc
real(r8) :: depth
real(r8) :: oerr

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities(routine)

! Read the namelist entry
call find_namelist_in_file("input.nml", "ascii_table_converter_nml", iunit)
read(iunit, nml = ascii_table_converter_nml, iostat = iocode)
call check_namelist_read(iunit, iocode, "ascii_table_converter_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=ascii_table_converter_nml)
if (do_nml_term()) write(     *     , nml=ascii_table_converter_nml)

! time setup
call set_calendar_type(GREGORIAN)
reference_time = set_date(2000,1,1) ! SMAP has times relative to midnight 1 Jan 2000
prev_time = set_time(0, 0)

! We need to know the maximum number of observations that might be added to
! the observation sequence.  The max possible number of obs needs to be 
! specified but it will only write out the actual number created.
! Each observation in this series will have a single
! observation value and a quality control flag.  
! Initialize two empty observations - one to track location
! in observation sequence - the other is for the new observation.

num_input_files = Check_Input_Files(input_file_list, filename_seq_list)

max_obs    = NOBS_PER_FILE*num_input_files
num_copies = 1
num_qc     = 1
first_obs  = .true.

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

if ( file_exist(obs_out_file) ) then ! existing file found, append to it

  call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)

else ! create a new one

  call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'observation')
  enddo
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  enddo
endif

!-----------------------------------------------------------------------
! Loop over all the input data files.

FileLoop: do ifile = 1,num_input_files

   filename = filename_seq_list(ifile)
   iunit    = open_file(filename, 'formatted', 'read')

   num_obs_read = 0
   num_obs_used = 0
   obsloop: do iline = 1,NOBS_PER_FILE    ! one observation per line
   
      read(iunit,*,iostat=iocode) timestamp, latitude, longitude, observation, soil_moisture_error, qc
      if (iocode < 0) exit obsloop
      if (iocode > 0) then
         write (string1,'(''Cannot read (error '',i3,'') line '',i8,'' in '',A)') &
                       iocode, iline, trim(filename)
         call error_handler(E_ERR,routine, string1, source)
      endif
      num_obs_read = num_obs_read + 1

      ! reject the observation for all kinds of reasons

      if (observation         <    0.0_r8) cycle obsloop  
!     if (soil_moisture_error <    0.0_r8) cycle obsloop  
      if (qc                  > maxgoodqc) cycle obsloop  

      ! convert the observation time to a DART time type

      days     = timestamp/86400    ! intentionally using integer arithmetic
      seconds  = timestamp - real(days,digits12)*86400.0_digits12
      offset   = set_time(seconds, days)
      obs_time = reference_time + offset
 
      if (num_obs_used < 1) then
         write(*,*)''
         write(*,*)'"'//trim(filename)//'"'
         write(*,*)'first good observation raw values: '
         call print_date(obs_time, ' date ')
         call print_time(obs_time, ' time ')
         write(*,*)'timestamp          ', timestamp
         write(*,*)'latitude           ', latitude
         write(*,*)'longitude          ', longitude
         write(*,*)'soil_moisture      ', observation
         write(*,*)'soil_moisture_error', soil_moisture_error
         write(*,*)'qc                 ', qc
   
         write(logfileunit,*)''
         write(logfileunit,*)'"'//trim(filename)//'"'
         write(logfileunit,*)'first good observation raw values: '
         call print_date(obs_time, ' date ', logfileunit)
         call print_time(obs_time, ' time ', logfileunit)
         write(logfileunit,*)'timestamp          ', timestamp
         write(logfileunit,*)'latitude           ', latitude
         write(logfileunit,*)'longitude          ', longitude
         write(logfileunit,*)'soil_moisture      ', observation
         write(logfileunit,*)'soil_moisture_error', soil_moisture_error
         write(logfileunit,*)'qc                 ', qc
      endif
   
      call get_time(obs_time, seconds, days)
    
      !>@todo put in a realistic observation error variance 
      oerr = soil_moisture_error   ! this one is all -99999ish at the present
      oerr = max(observation/10.0_r8, 0.02_r8)
      depth = 2.5_r8/100.0_r8  ! in meters

      if (longitude < 0.0_r8) longitude = longitude + 360.0_r8

      call create_3d_obs(latitude, longitude, depth, VERTISHEIGHT, observation, &
                         SOIL_MOISTURE, oerr, days, seconds, qc, obs)
      call add_obs_to_seq(obs_seq, obs, obs_time, prev_obs, prev_time, first_obs)
   
      num_obs_used = num_obs_used + 1

   enddo obsloop

   write(string1,*)'..  number of observations read ', num_obs_read
   write(string2,*)'number of good observations ', num_obs_used
   call error_handler(E_MSG, routine, string1, text2=string2)

enddo FileLoop   

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   write(string1,*)'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   if (verbose) call error_handler(E_MSG,routine,string1)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains

function Check_Input_Files(input_list, output_list)
character(len=*), intent(in)  :: input_list      !> filename containing list
character(len=*), intent(out) :: output_list(:)
integer                       :: Check_Input_Files

character(len=256) :: filename
character(len=256) :: ladjusted
integer :: iunit, iline, nlines

character(len=*), parameter :: routine='Check_Input_Files'

Check_Input_files = -1

call find_textfile_dims(input_list, nlines)

iunit = open_file(trim(input_list), 'formatted', 'read')

if (nlines >= MAX_NUM_INPUT_FILES ) then
   write(string1,*)'Too many files to process. Increase MAX_NUM_INPUT_FILES, recompile, and try again.'
   write(string2,*)'MAX_NUM_INPUT_FILES currently set to ',MAX_NUM_INPUT_FILES
   write(string3,*)'There were ',nlines,' files specified in ',trim(input_list)
   call error_handler(E_ERR,routine,string1,source, text2=string2, text3=string3)
endif

Check_Input_Files = 0
FileNameLoop: do iline = 1,nlines ! a lot of lines 

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=iocode) filename
   if (iocode > 0) then
      write(string1,*) 'While reading ', trim(input_list)
      write(string2,*) 'got read code (iostat) = ', iocode,' around line ',iline
      call error_handler(E_ERR, routine, string1, &
                    source, text2=string2)
   elseif (iocode < 0) then
      ! Normal end of file
      exit FileNameLoop
   else
      Check_Input_Files = Check_Input_Files + 1

      ladjusted = adjustl(filename)
      if ( file_exist(ladjusted) ) then
         output_list(Check_Input_Files) = ladjusted
      else
         write(string1,*)'following file does not exist:'
         call error_handler(E_ERR, routine, string1, &
             source, text2='"'//trim(ladjusted)//'"')
      endif
   endif

enddo FileNameLoop

end function Check_Input_Files

end program ascii_table_converter


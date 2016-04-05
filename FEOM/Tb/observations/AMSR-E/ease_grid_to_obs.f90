! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ease_grid_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   program to convert a series of flat binary files 
!
!   created 7 Nov 2013   Tim Hoar NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use          types_mod, only : r4, r8, PI, DEG2RAD

use      utilities_mod, only : initialize_utilities, finalize_utilities, &
                               open_file, close_file, find_namelist_in_file, &
                               check_namelist_read, nmlfileunit, get_unit, &
                               do_nml_file, do_nml_term, get_next_filename, &
                               error_handler, E_ERR, E_MSG, file_exist, &
                               find_textfile_dims

use   time_manager_mod, only : time_type, set_calendar_type, set_date, get_date, &
                               operator(>=), increment_time, set_time, get_time, &
                               operator(-), operator(+), GREGORIAN, &
                               print_time, print_date

use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                               static_init_obs_sequence, init_obs, write_obs_seq, & 
                               init_obs_sequence, get_num_obs, &
                               set_copy_meta_data, set_qc_meta_data

use            location_mod, only : VERTISSURFACE, set_location
use       obs_utilities_mod, only : add_obs_to_seq, create_3d_obs
use            obs_kind_mod, only : AMSRE_BRIGHTNESS_T
use obs_def_brightnessT_mod, only : set_amsre_metadata

use      EASE_utilities_mod, only : get_grid_dims, ezlh_inverse, read_ease_Tb, &
                                    deconstruct_filename, read_ease_TIM, &
                                    EASE_MISSING

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!----------------------------------------------------------------
! Namelist input with default values

character(len=256) :: input_file_list = 'file_list.txt'
character(len=256) :: obs_out_file = 'obs_seq.out'
logical            :: verbose = .false.
integer            :: max_num_input_files = 10

namelist /ease_grid_to_obs_nml/ &
         input_file_list, obs_out_file, verbose, max_num_input_files

!----------------------------------------------------------------

integer            :: num_input_files = 0  ! actual number of files
integer            :: ifile, istatus
character(len=256), allocatable, dimension(:) :: filename_seq_list
character(len=256) :: time_file_name

! information gleaned from filenaming convention
integer          :: iyear, idoy
character(len=2) :: gridarea
character(len=1) :: passdir
character(len=1) :: polarization
logical          :: is_time_file
real(r8)         :: footprint
real             :: frequency   ! real type to match EASE type

character(len=256) :: input_line
character(len=256) :: msgstring1,msgstring2,msgstring3

integer :: oday, osec, fday, fsec, iocode, iunit
integer :: num_copies, num_qc, max_obs
integer :: landcode = 0  ! FIXME ... totally bogus for now
           
logical  :: first_obs

! The EASE grid is a 2D array of 16 bit unsigned integers
integer, allocatable, dimension(:,:) :: Tb
integer, allocatable, dimension(:,:) :: Obs_Minutes
integer, dimension(2) :: ndims
integer :: key, nrows, ncols, irow, icol

real(r8) :: temp, terr, qc
real     :: rlat, rlon

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: file_time, time_obs, prev_time

!----------------------------------------------------------------
! start of executable code

call initialize_utilities('ease_grid_to_obs')

! time setup
call set_calendar_type(GREGORIAN)

! Read the namelist entry
call find_namelist_in_file("input.nml", "ease_grid_to_obs_nml", iunit)
read(iunit, nml = ease_grid_to_obs_nml, iostat = iocode)
call check_namelist_read(iunit, iocode, "ease_grid_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=ease_grid_to_obs_nml)
if (do_nml_term()) write(     *     , nml=ease_grid_to_obs_nml)

call find_textfile_dims(input_file_list, num_input_files)
if (num_input_files > max_num_input_files) then
   write(msgstring1,*) 'Found ',num_input_files,' input files in ',trim(input_file_list)
   write(msgstring2,*) 'Greater than "max_num_input_files" in input namelist (', max_num_input_files, ').'
   write(msgstring3,*) 'If you mean it, change "max_num_input_files" in the ease_grid_to_obs_nml.'
   call error_handler(E_ERR, 'main', msgstring1, source, revision, revdate, &
              text2=msgstring2, text3=msgstring3)
elseif (num_input_files > 0) then
   allocate(filename_seq_list(num_input_files))
else
   write(msgstring1,*) 'Found no valid filenames in ',trim(input_file_list)
   call error_handler(E_ERR, 'main', msgstring1, source, revision, revdate) 
endif

num_input_files = Check_Input_Files(input_file_list, filename_seq_list) 
write(*,*)' There are ',num_input_files,' input files.'

! need some basic information from the first file
istatus = deconstruct_filename( filename_seq_list(1), &
             gridarea, iyear, idoy, passdir, frequency, polarization, &
             is_time_file, time_file_name)
if (istatus /= 0) then
   write(msgstring2,*) 'filename nonconforming'
   call error_handler(E_ERR, 'main', trim(filename_seq_list(1)), &
                 source, revision, revdate, text2=msgstring2)
endif

ndims = get_grid_dims(gridarea)
nrows = ndims(1)
ncols = ndims(2)

allocate( Tb(nrows,ncols), Obs_Minutes(nrows,ncols) )

! each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
! There is a lot of observations per day. Should only do 1 day at at time.
max_obs    = nrows*ncols*num_input_files ! overkill
num_copies = 1
num_qc     = 1

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

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
call   set_qc_meta_data(obs_seq, 1,     'Data QC')

! --------------------------------------------------------------------------------
! Loop over all the input data files.

iunit = get_unit()  ! Get free unit for AMSR-E data file
FileLoop: do ifile = 1,num_input_files

   ! A little helpful logging
   write(msgstring1,*)'.. Converting file',ifile,' of ',num_input_files
   call error_handler(E_MSG, 'main', msgstring1, text2 = trim(filename_seq_list(ifile)))

   istatus = deconstruct_filename( filename_seq_list(ifile), &
             gridarea, iyear, idoy, passdir, frequency, polarization, &
             is_time_file, time_file_name)
   if (istatus /= 0) then
      call error_handler(E_ERR, 'main', 'filename nonconforming', &
             source, revision, revdate, text2= trim(filename_seq_list(ifile)))
   endif

   if (gridarea(2:2) == 'L') then
      footprint = 25.0_r8  ! EASE grid is 25km resolution
   else
      call error_handler(E_ERR, 'main', 'unknown footprint', &
             source, revision, revdate, text2= trim(filename_seq_list(ifile)))
   endif

   ! put date into a dart time format
   ! calculate first day of the year, then add in day-of-year from the file name
   file_time = set_date(iyear,1,1,0,0,0)
   call get_time(file_time, fsec, fday)
   file_time = set_time(fsec, fday+idoy-1)
   call get_time(file_time, fsec, fday)

   if (verbose) then
      write(*,*)trim(filename_seq_list(ifile))
      call print_date(file_time, str='file date is')
      call print_time(file_time, str='file time is')
   endif

   ! Read the time for each observation (minutes since filetime)
   iocode = read_ease_TIM(time_file_name, iunit, Obs_Minutes)

   ! Ignoring time files for now
   if ( is_time_file ) cycle FileLoop

   ! Fills up matrix of brightness temperatures
   iocode = read_ease_Tb(filename_seq_list(ifile), iunit, Tb)

   COLLOOP: do icol=1,ncols
   ROWLOOP: do irow=1,nrows

      if (          Tb(irow,icol) == 0            ) cycle ROWLOOP
      if ( Obs_Minutes(irow,icol) == EASE_MISSING ) cycle ROWLOOP

      ! Create the time of the observation by adding the base time
      ! for all the obs in the file to the number of minutes from the
      ! *.TIM file. Convert to DART time type and back to automatically
      ! make sure the number of seconds is [0,86400] etc.

      time_obs = set_time(fsec + Obs_Minutes(irow,icol)*60, fday)
      call get_time(time_obs, osec, oday)

      ! Convert icol,irow to lat/lon using EASE routine
      ! Intentional conversion between row-major and column-major nomenclature
      iocode = ezlh_inverse(gridarea, real(irow), real(icol), rlat, rlon)
      if (iocode /= 0) cycle ROWLOOP

      ! ensure the lat/lon values are in range
      if ( rlat >  90.0_r8 .or. rlat <  -90.0_r8 ) cycle ROWLOOP
      if ( rlon > 180.0_r8 .or. rlon < -180.0_r8 ) cycle ROWLOOP
      if ( rlon < 0.0_r8 ) rlon = rlon + 360.0_r8 ! changes into 0-360
   
      ! make an obs derived type, and then add it to the sequence
      temp = real(Tb(irow,icol),r8) / 10.0_r8
      terr = 2.0_r8   ! temps are [200,300] so this is about one percent

      call set_amsre_metadata(key, real(frequency,r8), footprint, polarization, landcode)

      call create_3d_obs(real(rlat,r8), real(rlon,r8), 0.0_r8, VERTISSURFACE, temp, &
                            AMSRE_BRIGHTNESS_T, terr, oday, osec, qc, obs, key)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
   enddo ROWLOOP
   enddo COLLOOP

enddo FileLoop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

deallocate(Tb,filename_seq_list)

! end of main program
call finalize_utilities()


contains


function Check_Input_Files(input_list, output_list)
! Read a list of files to process 
character(len=*),               intent(in)  :: input_list
character(len=*), dimension(:), intent(out) :: output_list
integer                                     :: Check_Input_Files

character(len=256) :: ladjusted
integer :: iline, fnamelen

Check_Input_files = -1

iunit = open_file(trim(input_list), 'formatted', 'read')

Check_Input_Files = 0
FileNameLoop: do iline = 1,size(output_list) ! a lot of lines 

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=iocode) input_line
   if (iocode > 0) then 
      write(msgstring1,*) 'While reading ', trim(input_list)
      write(msgstring2,*) 'got read code (iostat) = ', iocode,' around line ',iline
      call error_handler(E_ERR, 'Check_Input_Files', msgstring1, &
                    source, revision, revdate, text2=msgstring2)
   elseif (iocode < 0) then 
      ! Normal end of file
      exit FileNameLoop
   else

      ladjusted = adjustl(input_line)
      fnamelen  = len_trim(ladjusted)

      ! Remove any file with a .TIM extension.
      ! These are TIME files, not data files.

      if ( ladjusted(fnamelen-2:fnamelen) == 'TIM' ) cycle FileNameLoop

      Check_Input_Files = Check_Input_Files + 1

      if ( file_exist(trim(ladjusted)) ) then
         output_list(Check_Input_Files) = trim(ladjusted)
      else
         write(msgstring1,*)'file does not exist.'
         call error_handler(E_ERR,'Check_Input_Files',&
          msgstring1,source,revision,revdate,text2=trim(ladjusted))
      endif
   endif

enddo FileNameLoop

end function Check_Input_Files


end program ease_grid_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

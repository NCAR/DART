! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ease_grid_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   program to convert a series of HDF5 files
!
!   created 13 Nov 2017   Tim Hoar NCAR/IMAGe
!
! https://nsidc.org/data/ease/tools
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! /glade/u/home/shimj/DART/observations/obs_converters/EASE-Grid/SMAP_retrival.log

use          types_mod, only : r4, r8, digits12

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

use       location_mod, only : VERTISHEIGHT, set_location

use  obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use  netcdf_utilities_mod, only : nc_get_variable, nc_check

use       obs_kind_mod, only : SOIL_MOISTURE

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

!----------------------------------------------------------------
! Namelist input with default values

character(len=256) :: input_file_list = 'file_list.txt'
character(len=256) :: obs_out_file = 'obs_seq.out'
logical            :: verbose = .false.

namelist /ease_grid_to_obs_nml/ &
         input_file_list, obs_out_file, verbose

!----------------------------------------------------------------

! MAX_NUM_INPUT_FILES : max number of input files to be processed
integer, parameter :: MAX_NUM_INPUT_FILES = 500
integer            :: num_input_files = 0  ! actual number of files
integer            :: ifile, istatus
character(len=256), dimension(MAX_NUM_INPUT_FILES) :: filename_seq_list
character(len=256) :: filename

! information gleaned from filenaming convention
integer          :: iyear, idoy
character(len=2) :: gridarea
character(len=1) :: passdir
character(len=1) :: polarization
logical          :: is_time_file
real(r8)         :: footprint
real             :: frequency   ! real type to match EASE type

character(len=256) :: input_line
character(len=512) :: string1, string2, string3

integer :: oday, osec, iocode, iunit, otype
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs
integer :: landcode = 0  ! FIXME ... totally bogus for now
           
logical  :: first_obs

! The EASE grid 
real(r8),        allocatable, dimension(:) :: longitude
real(r8),        allocatable, dimension(:) :: latitude
real(r8),        allocatable, dimension(:) :: observation
real(r8),        allocatable, dimension(:) :: soil_moisture_error
type(time_type), allocatable, dimension(:) :: obs_time
integer,         allocatable, dimension(:) :: retrieval_flag

integer                      :: dimids(NF90_MAX_DIMS)
integer                      :: dimlens(NF90_MAX_DIMS)
character(len=NF90_MAX_NAME) :: dimnames(NF90_MAX_DIMS)

integer :: key, icount, ncid, VarID, io
integer :: counts = 20000

character(len=NF90_MAX_NAME) :: varname
integer  :: xtype, ndims, nAtts

real(r8) :: temp, terr, qc
real(r8) :: rlat, rlon, depth_cm, depth_m

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: cal_day0, observation_time, prev_time

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

num_input_files = Check_Input_Files(input_file_list, filename_seq_list) 

! each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
! There is a lot of observations per day. Should only do 1 day at at time.
max_obs    = counts*num_input_files ! overkill
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

FileLoop: do ifile = 1,num_input_files

   filename = filename_seq_list(ifile)

   ! A little helpful logging
   write(string1,*)'.. Converting file',ifile,' of ',num_input_files
   call error_handler(E_MSG, 'ease_grid_to_dart', string1, text2 = trim(filename))

   io = nf90_open(trim(filename), nf90_nowrite, ncid)
   call nc_check(io, 'ease_grid_to_dart', context='nf90_open',filename=filename)

   ! Get dimension information from the longitude variable

   io = nf90_inq_varid(ncid,'longitude',VarID)
   call nc_check(io, 'ease_grid_to_dart', context='inq_varid longitude', filename=filename)

   io = nf90_inquire_variable(ncid, VarID, varname, xtype, ndims, dimids, nAtts)
   call nc_check(io, 'ease_grid_to_dart', context='inquire_variable longitude', filename=filename)

   !>@todo FIXME check to make sure there is only 1 dimension

   io = nf90_inquire_dimension(ncid, dimids(1), len=dimlens(1))
   call nc_check(io, 'ease_grid_to_dart', context='inquire_dimension longitude', filename=filename)

   counts = dimlens(1) 

   allocate(          longitude(counts))
   allocate(           latitude(counts))
   allocate(      observation(counts)) 
   allocate(soil_moisture_error(counts))  
   allocate(           obs_time(counts))  
   allocate(     retrieval_flag(counts))  

   call nc_get_variable(ncid, 'longitude', longitude, 'ease_grid_to_dart', filename)
   call nc_get_variable(ncid, 'latitude', latitude, 'ease_grid_to_dart', filename)
   call nc_get_variable(ncid, 'soil_moisture', observation, 'ease_grid_to_dart', filename)
   call nc_get_variable(ncid, 'soil_moisture_error', soil_moisture_error, 'ease_grid_to_dart', filename)
   call nc_get_variable(ncid, 'retrieval_qual_flag', retrieval_flag, 'ease_grid_to_dart', filename)
   call read_observation_times(ncid, filename)

   COUNTLOOP: do icount=1,counts

      if ( observation(icount) ==  0.0_r8 ) cycle COUNTLOOP

      rlat = latitude(icount)
      rlon = longitude(icount)

      ! ensure the lat/longitude values are in range
      if ( rlat >  90.0_r8 .or. rlat <  -90.0_r8 ) cycle COUNTLOOP
      if ( rlon > 360.0_r8 .or. rlon <    0.0_r8 ) cycle COUNTLOOP

      depth_cm = 5
      depth_m = real(depth_cm,r8)/100.0_r8

      call get_time(obs_time(icount), osec, oday)
      qc = retrieval_flag(icount)

      call create_3d_obs(rlat, rlon, depth_m/1000.0_r8, VERTISHEIGHT, observation(icount), &
                        SOIL_MOISTURE, soil_moisture_error(icount), oday, osec, qc, obs)

      call add_obs_to_seq(obs_seq, obs, obs_time(icount), prev_obs, prev_time, first_obs)

   enddo COUNTLOOP

   deallocate(longitude)
   deallocate(latitude)
   deallocate(observation)
   deallocate(soil_moisture_error)
   deallocate(obs_time)
   deallocate(retrieval_flag)

enddo FileLoop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains

!-----------------------------------------------------------------------
!> Read a list of files to process.
!> Make sure each of the files exists.

function Check_Input_Files(input_list, output_list)
character(len=*), intent(in)  :: input_list      !> filename containing list
character(len=*), intent(out) :: output_list(:)
integer                       :: Check_Input_Files

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
   call error_handler(E_ERR,routine,string1,source,revision,revdate, text2=string2, text3=string3)
endif

Check_Input_Files = 0
FileNameLoop: do iline = 1,nlines ! a lot of lines 

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=iocode) input_line
   if (iocode > 0) then 
      write(string1,*) 'While reading ', trim(input_list)
      write(string2,*) 'got read code (iostat) = ', iocode,' around line ',iline
      call error_handler(E_ERR, routine, string1, &
                    source, revision, revdate, text2=string2)
   elseif (iocode < 0) then 
      ! Normal end of file
      exit FileNameLoop
   else
      Check_Input_Files = Check_Input_Files + 1

      ladjusted = adjustl(input_line)
      if ( file_exist(trim(ladjusted)) ) then
         output_list(Check_Input_Files) = trim(ladjusted)
      else
         write(string1,*)'file does not exist.'
         call error_handler(E_ERR,routine,&
          string1,source,revision,revdate,text2=ladjusted)
      endif
   endif

enddo FileNameLoop

end function Check_Input_Files



subroutine read_observation_times(ncid,filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename

real(digits12) :: seconds_array(size(obs_time))
real(digits12) :: remainder
type(time_type) :: base_time, offset
integer :: days, seconds, itime

character(len=512) :: units
character(len=512) :: long_name

character(len=*), parameter :: routine = 'read_observation_times'

call nc_get_variable(ncid, 'tb_time_seconds', seconds_array, routine,filename )

io = nf90_inq_varid(ncid, 'tb_time_seconds', VarID)
call nc_check(io, routine, 'inq_varid tb_time_seconds', filename=filename)

io = nf90_get_att(ncid, VarID, 'units', units)
call nc_check(io, routine, 'get_att tb_time_seconds units', filename=filename)

write(*,*)'TJH units is ',trim(units)

io = nf90_get_att(ncid, VarID, 'long_name', long_name)
call nc_check(io, routine, 'get_att tb_time_seconds long_name', filename=filename)

write(*,*)'TJH long_name is ',trim(long_name)

! hardcoded for now ... could check long_name
base_time = set_date(2000,1,1,0,0,0)

do itime = 1,size(obs_time)
   days = seconds_array(itime)/86400
   remainder = seconds_array(itime) - real(days,digits12) * 86400.0_digits12
   seconds = floor(remainder*86400.0_digits12)
   offset = set_time(seconds,days)
   obs_time(itime) = base_time + offset
enddo

if (verbose) then
   call print_date(obs_time(1), str='date is')
   call print_time(obs_time(1), str='time is')
endif

end subroutine read_observation_times

end program ease_grid_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

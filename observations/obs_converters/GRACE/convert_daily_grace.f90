! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program convert_daily_grace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convert_daily_grace - reads daily GRACE TOTAL WATER STORAGE 
!                       observations and writes a DART obs_seq file.
!
! One day in, one day out
!
! The input data to this routine have been preprocessed to convert
! the monthly GRACE total water storage anomalies to DAILY total water storage.
! Part of the preprocessing is to convert the monthly anomalies to daily
! anomalies, and then it is required to add in the long-term mean total water
! storage from whatever model you are using.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8

use      location_mod, only : VERTISUNDEF

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              nmlfileunit, do_nml_file, do_nml_term, &
                              find_namelist_in_file, check_namelist_read, &
                              error_handler, E_ERR, E_MSG, find_textfile_dims, &
                              open_file, close_file, file_exist

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, print_date, get_time, print_time

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data, &
                              destroy_obs_sequence

use      obs_kind_mod, only : GRACE_TOTAL_WATER_STORAGE

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use netcdf_utilities_mod, only : nc_check, nc_open_file_readonly, nc_close_file, &
                                 nc_get_variable_size, nc_get_variable

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"
character(len=*), parameter :: routine  = 'convert_daily_grace'

character(len=512) :: string1,string2,string3

!----------------------------------------------------------------
! Namelist input with default values

character(len=256) :: input_file_list = 'file_list.txt'
logical            :: verbose = .false.

namelist /convert_daily_grace_nml/ input_file_list, verbose

!----------------------------------------------------------------

integer, parameter :: MAX_NUM_INPUT_FILES = 500
integer            :: num_input_files = 0
integer            :: ifile

character(len=256), dimension(MAX_NUM_INPUT_FILES) :: filename_seq_list

character(len=256) :: output_name

integer :: iocode, iunit
integer :: osec, oday
integer :: year, month, day
integer :: num_copies, num_qc
           
logical  :: first_obs

real(r8), allocatable, dimension(:) :: TWS
real(r8), allocatable, dimension(:) :: ERR
real(r8), allocatable, dimension(:) :: LON
real(r8), allocatable, dimension(:) :: LAT

integer :: icount, counts, ncid

real(r8) :: qc
real(r8) :: rlat, rlon

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

!----------------------------------------------------------------
! start of executable code

call initialize_utilities(routine)

! time setup
call set_calendar_type(GREGORIAN)

! Read the namelist entry
call find_namelist_in_file("input.nml", "convert_daily_grace_nml", iunit)
read(iunit, nml = convert_daily_grace_nml, iostat = iocode)
call check_namelist_read(iunit, iocode, "convert_daily_grace_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_daily_grace_nml)
if (do_nml_term()) write(     *     , nml=convert_daily_grace_nml)

num_input_files = Check_Input_Files(input_file_list, filename_seq_list) 
write(*,*)' There are ',num_input_files,' input files.'

! call the initialization code, and initialize two empty observation types
! each observation in this series will have a single observation value 
! and a quality control flag.  Set the DART data quality control, 0 is good. 
! increasingly larger QC values are more questionable quality data

num_copies = 1
num_qc     = 1
qc         = 0.0_r8

call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! --------------------------------------------------------------------------------
! Loop over all the input data files.

FileLoop: do ifile = 1,num_input_files
   
   first_obs = .true.

   ! Fills up matrix of grace tws
   ncid = nc_open_file_readonly(filename_seq_list(ifile), routine)
   
   call nc_get_variable_size(ncid, 'tws', counts)

   allocate( TWS(counts) )
   allocate( ERR(counts) )
   allocate( LON(counts))
   allocate( LAT(counts))

   call nc_get_variable(ncid, 'lat', LAT, routine)
   call nc_get_variable(ncid, 'lon', LON, routine)
   call nc_get_variable(ncid, 'tws', TWS, routine)
   call nc_get_variable(ncid, 'err', ERR, routine)

   where(LON < 0.0_r8) LON = LON + 360.0_r8

   call nc_close_file(ncid, routine)
   
   ! create a new, empty obs_seq file.  you must give a max limit
   ! on number of obs.  increase the size if too small.
   call init_obs_sequence(obs_seq, num_copies, num_qc, counts)

   ! the first one needs to contain the string 'observation' and the
   ! second needs the string 'QC'.
   call set_copy_meta_data(obs_seq, 1, 'observation')
   call   set_qc_meta_data(obs_seq, 1,     'Data QC')
   
   call get_yyyymmdd(filename_seq_list(ifile), year, month, day)

   write(output_name,'(''obs_seq.'',i4.4,''-'',i2.2,''-'',i2.2,''.out'')') &
                     year, month, day
 
   time_obs = set_date(year,month,day,0,0,0)
   call get_time(time_obs, osec, oday)

   ! A little helpful logging
   if (verbose) then
      write(*,*)'converting "',trim(filename_seq_list(ifile)),'"'
      write(string1,*)'"'//trim(output_name)//'" date is'
      write(string2,*)'"'//trim(output_name)//'" time is'
      call print_date(time_obs, str=string1)
      call print_time(time_obs, str=string2)
   endif

   COUNTLOOP: do icount=1,counts

      if ( TWS(icount) <  0.0_r8 ) cycle COUNTLOOP

      rlat = LAT(icount)
      rlon = LON(icount)

      ! ensure the lat/lon values are in range
      if ( rlat >  90.0_r8 .or. rlat <  -90.0_r8 ) cycle COUNTLOOP
      if ( rlon > 360.0_r8 .or. rlon <    0.0_r8 ) cycle COUNTLOOP

      ! make an obs derived type, and then add it to the sequence
      
      call create_3d_obs(rlat, rlon, 0.0_r8, VERTISUNDEF, TWS(icount), &
                  GRACE_TOTAL_WATER_STORAGE, ERR(icount), oday, osec, qc, obs)

      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

   enddo COUNTLOOP


   ! if we added any obs to the sequence, write it out to a file now.
   if ( get_num_obs(obs_seq) > 0 ) then
      if (verbose) then
         write(string1,*)'"'//trim(output_name)//'" has ', &
                        get_num_obs(obs_seq), ' observations.'
         call error_handler(E_MSG, routine, string1)
      endif
      call write_obs_seq(obs_seq, output_name)
   endif

   deallocate(TWS)
   deallocate(ERR)
   deallocate(LON)
   deallocate(LAT)

   call destroy_obs_sequence(obs_seq)

enddo FileLoop

! end of main program
call finalize_utilities()

contains


function Check_Input_Files(input_list, output_list)
! Read a list of files to process 
character(len=*),               intent(in)  :: input_list
character(len=*), dimension(:), intent(out) :: output_list
integer                                     :: Check_Input_Files

character(len=256) :: input_line
character(len=256) :: ladjusted
integer, parameter :: MAXLINES = 500             ! Long; change from 1000
integer :: iline

iunit = open_file(trim(input_list), 'formatted', 'read')

Check_Input_Files = 0
FileNameLoop: do iline = 1,MAXLINES ! a lot of lines 

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=iocode) input_line
   if (iocode > 0) then 
      write(string1,*) 'While reading ', trim(input_list)
      write(string2,*) 'got read code (iostat) = ', iocode,' around line ',iline
      call error_handler(E_ERR, 'Check_Input_Files', string1, &
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
         call error_handler(E_ERR,'Check_Input_Files',&
          string1,source,revision,revdate,text2=trim(ladjusted))
      endif
   endif

enddo FileNameLoop

if (Check_Input_Files >= MAXLINES-1 ) then
   write(string1,*)'Too many files to process. Increase MAXLINES and try again.'
   call error_handler(E_ERR,'Check_Input_Files',string1,source,revision,revdate)
endif

end function Check_Input_Files


!> get the year,month,day from the name of input files

subroutine get_yyyymmdd(filename, year, month, day)
character(len=*), intent(in) :: filename
integer,          intent(out) :: year, month, day

! The filenames END with "GRACE_1d_2003-01-01"
! There is no file extension.

integer :: strlen, index1, io

strlen = len_trim(filename)
index1 = strlen - 10 + 1  !    YYYY-MM-DD is 10 chars

read(filename(index1:strlen),"(I4,1X,I2,1X,I2)",iostat=io) year, month, day   
if (io /= 0) then
   write(string1,*)'Unable to read year, month, day from filename.'
   write(string2,*)'filename was "'//trim(filename)//'"'
   write(string3,*)'Expected last 10 chars to contain the date, they are "' &
                   & //filename(index1:strlen)//'"'
   call error_handler(E_ERR, 'get_yyymmdd', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

end subroutine get_yyyymmdd
  
 
end program convert_daily_grace 

! <next few lines under version control, do not edit>
! $URL$
! $Revision$
! $Date$


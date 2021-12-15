! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program SMAP_L2_to_obs

!=======================================================================
!  program to convert a series of HDF5 files
!
!   created 13 Nov 2017   Tim Hoar NCAR/IMAGe
!
! https://nsidc.org/data/ease/tools
! https://nsidc.org/data/SPL2SMP/versions/4
!
! "This Level-2 (L2) soil moisture product provides estimates of global land 
! surface conditions retrieved by the Soil Moisture Active Passive (SMAP) 
! passive microwave radiometer during 6:00 a.m. descending and 6:00 p.m. 
! ascending half-orbit passes. SMAP L-band brightness temperatures are used  to
! derive soil moisture data, which are then resampled to an Earth-fixed, global,
! cylindrical 36 km Equal-Area Scalable Earth Grid, Version 2.0 (EASE-Grid 2.0).
!
! Data Set ID: SPL2SMP
! SMAP L2 Radiometer Half-Orbit 36 km EASE-Grid Soil Moisture, Version 4
!
! "Surface soil moisture (0-5 cm) in m3/m3 derived from brightness temperatures
! (TBs) is output on a fixed global 36 km EASE-Grid 2.0. Also included are 
! brightness temperatures in kelvin representing the weighted average of 
! Level-1B brightness temperatures whose boresights fall within a 36 km 
! EASE-Grid 2.0 cell."
!
!=======================================================================

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

use HDF5_utilities_mod, only : h5_open, h5_check, &
                               h5_get_rank, &
                               h5_get_dimensions, &
                               h5_get_dset_dspace

use       obs_kind_mod, only : SOIL_MOISTURE

use HDF5
use H5LT

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'SMAP_L2_to_obs.f90'

character(len=*), parameter :: routine = 'SMAP_L2_to_obs'

!-----------------------------------------------------------------------
! Namelist input with default values

character(len=256) :: input_file_list = 'file_list.txt'
character(len=256) :: obs_out_file = 'obs_seq.out'
logical            :: verbose = .false.

namelist /SMAP_L2_to_obs_nml/ &
         input_file_list, obs_out_file, verbose

!-----------------------------------------------------------------------

! MAX_NUM_INPUT_FILES : max number of input files to be processed
integer, parameter :: MAX_NUM_INPUT_FILES = 500
integer            :: num_input_files = 0  ! actual number of files
integer            :: ifile
character(len=256), dimension(MAX_NUM_INPUT_FILES) :: filename_seq_list
character(len=256) :: filename

character(len=512) :: string1, string2, string3

integer :: oday, osec, iocode, iunit
integer :: num_copies, num_qc, max_obs
           
logical :: first_obs

! The EASE grid 
real(r4),        allocatable, dimension(:) :: longitude
real(r4),        allocatable, dimension(:) :: latitude
real(r4),        allocatable, dimension(:) :: observation
real(r4),        allocatable, dimension(:) :: soil_moisture_error_std
type(time_type), allocatable, dimension(:) :: obs_time
integer,         allocatable, dimension(:) :: retrieval_flag

integer :: icount
integer :: counts = 20000

character(len=256) :: varname ! HDF variable names can be long
integer  :: ndims

real(r8) :: qc, obs_val, err_std
real(r8) :: rlat, rlon, depth_cm, depth_m
real(r4) :: obs_FillValue(1), obs_valid_min(1), obs_valid_max(1)
real(r4) :: sme_FillValue(1), sme_valid_min(1), sme_valid_max(1)
real(r4) :: rqf_FillValue(1)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time

integer(HSIZE_T), allocatable :: dimlens(:)
integer(HID_T) :: file_id, dset_id, dspace_id
integer        :: hdferr

!-----------------------------------------------------------------------
! start of executable code

call initialize_utilities(routine)

! time setup
call set_calendar_type(GREGORIAN)

! Read the namelist entry
call find_namelist_in_file("input.nml", "SMAP_L2_to_obs_nml", iunit)
read(iunit, nml = SMAP_L2_to_obs_nml, iostat = iocode)
call check_namelist_read(iunit, iocode, "SMAP_L2_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=SMAP_L2_to_obs_nml)
if (do_nml_term()) write(     *     , nml=SMAP_L2_to_obs_nml)

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

!-----------------------------------------------------------------------
! Loop over all the input data files.

FileLoop: do ifile = 1,num_input_files

   filename = filename_seq_list(ifile)

   ! A little helpful logging
   write(string1,*)'.. Converting file',ifile,' of ',num_input_files
   call error_handler(E_MSG, routine, string1, text2 = '"'//trim(filename)//'"')

   file_id = h5_open(filename, H5F_ACC_RDONLY_F)

   ! Get dimension information from the longitude variable

   varname = '/Soil_Moisture_Retrieval_Data/longitude'

   write(string1,*) routine, trim(filename), trim(varname)

   call h5_get_dset_dspace(file_id, varname, dset_id, dspace_id, string1)
   ndims = h5_get_rank(dspace_id, string1)

   if (ndims /= 1) then
      call error_handler(E_ERR,routine,'wrong shape', source)
   endif  

   allocate(dimlens(ndims)) 
   call h5_get_dimensions(dspace_id, dimlens, context=string1)

   counts = dimlens(1) 

   allocate(          longitude(counts))
   allocate(           latitude(counts))
   allocate(        observation(counts)) 
   allocate(soil_moisture_error_std(counts))  
   allocate(           obs_time(counts))  
   allocate(     retrieval_flag(counts))  

   call read_observation_times(file_id,filename)

   varname = '/Soil_Moisture_Retrieval_Data/longitude'
   call h5ltread_dataset_float_f(file_id,varname,longitude,dimlens,hdferr)
   call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

   varname = '/Soil_Moisture_Retrieval_Data/latitude'
   call h5ltread_dataset_float_f(file_id,varname,latitude,dimlens,hdferr)
   call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

   varname = '/Soil_Moisture_Retrieval_Data/soil_moisture'
   call h5ltread_dataset_float_f(file_id,varname,observation,dimlens,hdferr)
   call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'_FillValue',obs_FillValue,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'valid_min',obs_valid_min,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'valid_max',obs_valid_max,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   !>@todo The soil moisture error values are ALL MISSING in the data files I have.
   !> not using this for now.
   varname = '/Soil_Moisture_Retrieval_Data/soil_moisture_error'
   call h5ltread_dataset_float_f(file_id,varname,soil_moisture_error_std,dimlens,hdferr)
   call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'_FillValue',sme_FillValue,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'valid_min',sme_valid_min,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'valid_max',sme_valid_max,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   varname = '/Soil_Moisture_Retrieval_Data/retrieval_qual_flag'
   call h5ltread_dataset_int_f(file_id,varname,retrieval_flag,dimlens,hdferr)
   call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'_FillValue',rqf_FillValue,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   if (verbose) then
      write(*,*)
      write(*,*)' obs FillValue is ',obs_FillValue
      write(*,*)' obs valid_min is ',obs_valid_min
      write(*,*)' obs valid_max is ',obs_valid_max
      write(*,*)' sme FillValue is ',sme_FillValue
      write(*,*)' sme valid_min is ',sme_valid_min
      write(*,*)' sme valid_max is ',sme_valid_max
      write(*,*)' rqf FillValue is ',rqf_FillValue
      write(*,*)
   endif

   COUNTLOOP: do icount=1,counts

      if (            observation(icount) == obs_FillValue(1)) cycle COUNTLOOP
!     if (soil_moisture_error_std(icount) == sme_FillValue(1)) cycle COUNTLOOP
      if (         retrieval_flag(icount) == rqf_FillValue(1)) cycle COUNTLOOP

      if (            observation(icount) < obs_valid_min(1)) cycle COUNTLOOP
      if (            observation(icount) > obs_valid_max(1)) cycle COUNTLOOP
!     if (soil_moisture_error_std(icount) < sme_valid_min(1)) cycle COUNTLOOP
!     if (soil_moisture_error_std(icount) > sme_valid_max(1)) cycle COUNTLOOP

      rlat = real(latitude( icount),r8)
      rlon = real(longitude(icount),r8)

      ! ensure the lat/longitude values are in range
      if ( rlat >  90.0_r8 .or. rlat <  -90.0_r8 ) cycle COUNTLOOP
      if ( rlon > 360.0_r8 .or. rlon <    0.0_r8 ) cycle COUNTLOOP

      depth_cm = 5.0_r8/2.0_r8
      depth_m = real(depth_cm,r8)/100.0_r8

      call get_time(obs_time(icount), osec, oday)
      qc = retrieval_flag(icount)

      obs_val = real(observation(icount),r8)

      ! Since the soil moisture error values are all _FillValue in the data
      ! I have, we are simply assuming something like 20% ... 
      ! cannot have an error standard deviation of 0.0, so there must be an alternative
      ! minimum error specification ... just using 0.01 until proven otherwise
      ! Both 20% and 0.01 have no scientific basis and should be explored.

      ! err_std = real(soil_moisture_error_std(icount),r8)

      err_std = max(obs_val/20.0_r8, 0.01_r8)
 
      call create_3d_obs(rlat, rlon, depth_m, VERTISHEIGHT, obs_val, &
                        SOIL_MOISTURE, err_std, oday, osec, qc, obs)

      call add_obs_to_seq(obs_seq, obs, obs_time(icount), prev_obs, prev_time, first_obs)

   enddo COUNTLOOP

   deallocate(dimlens)
   deallocate(longitude)
   deallocate(latitude)
   deallocate(observation)
   deallocate(soil_moisture_error_std)
   deallocate(obs_time)
   deallocate(retrieval_flag)

enddo FileLoop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   write(string1,*) 'writing observations: obs_count = ', get_num_obs(obs_seq)
   call error_handler(E_MSG, routine, string1)
   call write_obs_seq(obs_seq, obs_out_file)
else
   call error_handler(E_MSG, routine, 'no observations to write out.')
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
   call error_handler(E_ERR,routine,string1,source,text2=string2, text3=string3)
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
      if ( file_exist(trim(ladjusted)) ) then
         output_list(Check_Input_Files) = ladjusted
      else
         write(string1,*)'following file does not exist:'
         call error_handler(E_ERR, routine, string1, &
             source, text2='"'//trim(ladjusted)//'"')
      endif
   endif

enddo FileNameLoop

end function Check_Input_Files


!-----------------------------------------------------------------------
!> read the observation time array and convert to an array of DART
!> time_type objects

subroutine read_observation_times(file_id, filename)

integer(HID_T),   intent(in) :: file_id
character(len=*), intent(in) :: filename


real(digits12) :: seconds_array(size(obs_time))
real(digits12) :: remainder
type(time_type) :: base_time, offset
integer :: days, seconds, itime

character(len=512) :: units
character(len=512) :: long_name

character(len=*), parameter :: routine = 'read_observation_times'

varname = '/Soil_Moisture_Retrieval_Data/tb_time_seconds'
call h5ltread_dataset_double_f(file_id, varname, seconds_array, dimlens, hdferr)
call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

call h5ltget_attribute_string_f(file_id, varname, 'units', units, hdferr)
call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

!>@todo  check to make sure units are 'seconds'
! write(*,*)'TJH units is ',trim(units)

call h5ltget_attribute_string_f(file_id, varname, 'long_name', long_name, hdferr)
call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

!>@todo  check to make sure 'seconds since midnight on January 1, 2000 UTC' 
! is somewhere in the long_name
! write(*,*)'TJH long_name is ',trim(long_name)

!>@todo hardcoded for now 
base_time = set_date(2000,1,1,0,0,0)

do itime = 1,size(obs_time)
   days      = seconds_array(itime)/86400
   remainder = seconds_array(itime) - real(days,digits12) * 86400.0_digits12
   seconds   = floor(remainder)
   offset    = set_time(seconds,days)
   obs_time(itime) = base_time + offset
enddo

if (verbose) then
   call print_date(obs_time(1), str='date is')
   call print_time(obs_time(1), str='time is')
endif

end subroutine read_observation_times

end program SMAP_L2_to_obs


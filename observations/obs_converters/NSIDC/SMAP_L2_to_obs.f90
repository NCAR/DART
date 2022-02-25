! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
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

use       obs_kind_mod, only : SOIL_MOISTURE, SOIL_TEMPERATURE

use HDF5
use H5LT

implicit none

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
logical :: moisture_good, temperature_good

! The EASE grid
real(r4),        allocatable, dimension(:) :: longitude
real(r4),        allocatable, dimension(:) :: latitude
real(r4),        allocatable, dimension(:) :: moisture
real(r4),        allocatable, dimension(:) :: moisture_error_std
real(r4),        allocatable, dimension(:) :: temperature
type(time_type), allocatable, dimension(:) :: obs_time
integer,         allocatable, dimension(:) :: retrieval_flag

integer :: icount
integer :: counts = 40000

integer :: moisture_fill_count
integer :: moisture_min_count
integer :: moisture_max_count
integer :: moisture_count
integer :: temperature_fill_count
integer :: temperature_min_count
integer :: temperature_max_count
integer :: temperature_count
integer :: retrieval_flag_count
integer :: badlat_count
integer :: badlon_count

character(len=256) :: varname ! HDF variable names can be long
integer  :: ndims

real(r8) :: qc, obs_val, err_std
real(r8) :: rlat, rlon, depth_cm, sm_depth_m, t_depth_m
real(r4) :: temperature_FillValue(1), temperature_valid_min(1), temperature_valid_max(1)
real(r4) ::    moisture_FillValue(1),    moisture_valid_min(1),    moisture_valid_max(1)
real(r4) ::         sme_FillValue(1),         sme_valid_min(1),         sme_valid_max(1)
real(r4) ::         rqf_FillValue(1)

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

! Purely empirical stab at the depth of the soil moisture observation
depth_cm = 5.0_r8/2.0_r8
sm_depth_m = real(depth_cm,r8)/100.0_r8

! Purely empirical stab at the depth of the 'surface_temperature'
depth_cm = 1.0_r8/2.0_r8
t_depth_m = real(depth_cm,r8)/100.0_r8

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

   allocate(         longitude(counts))
   allocate(          latitude(counts))
   allocate(          moisture(counts))
   allocate(       temperature(counts))
   allocate(moisture_error_std(counts))
   allocate(          obs_time(counts))
   allocate(    retrieval_flag(counts))

   call read_observation_times(file_id,filename)

   varname = '/Soil_Moisture_Retrieval_Data/longitude'
   call h5ltread_dataset_float_f(file_id,varname,longitude,dimlens,hdferr)
   call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

   varname = '/Soil_Moisture_Retrieval_Data/latitude'
   call h5ltread_dataset_float_f(file_id,varname,latitude,dimlens,hdferr)
   call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

   varname = '/Soil_Moisture_Retrieval_Data/soil_moisture'
   call h5ltread_dataset_float_f(file_id,varname,moisture,dimlens,hdferr)
   call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'_FillValue',moisture_FillValue,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'valid_min',moisture_valid_min,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'valid_max',moisture_valid_max,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   ! From /gpfs/fs1/work/thoar/DART/cesm_clm/observations/obs_converters/NSIDC
   ! 'It has come to our attention that this parameter is often mistaken for the 
   !  physical temperature of the top soil layer. The designation "effective" signifies
   !  an attempt to capture the soil integrated temperature and canopy temperature in a
   !  single parameter, as is widely reported in the literature.  Depending on the actual
   !  emission sensing depth (which varies with soil moisture), this parameter usually
   !  does not coincide with a thermal physical temperature at a fixed depth
   ! (e.g. 5 cm or 10 cm).' bummer ...

   varname = '/Soil_Moisture_Retrieval_Data/surface_temperature'
   call h5ltread_dataset_float_f(file_id,varname,temperature,dimlens,hdferr)
   call h5_check(hdferr,routine,'h5ltread_dataset_float_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'_FillValue',temperature_FillValue,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'valid_min',temperature_valid_min,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   call h5ltget_attribute_float_f(file_id,varname,'valid_max',temperature_valid_max,hdferr)
   call h5_check(hdferr,routine,'h5ltget_attribute_string_f',varname,filename)

   !>@todo The soil moisture error values are ALL MISSING in the data files I have.
   !> not using this for now.
   varname = '/Soil_Moisture_Retrieval_Data/soil_moisture_error'
   call h5ltread_dataset_float_f(file_id,varname,moisture_error_std,dimlens,hdferr)
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
      write(*,*)' moisture  FillValue is ',moisture_FillValue
      write(*,*)' moisture  valid_min is ',moisture_valid_min
      write(*,*)' moisture  valid_max is ',moisture_valid_max
      write(*,*)' surf_temp FillValue is ',temperature_FillValue
      write(*,*)' surf_temp valid_min is ',temperature_valid_min
      write(*,*)' surf_temp valid_max is ',temperature_valid_max
      write(*,*)' sme       FillValue is ',sme_FillValue
      write(*,*)' sme       valid_min is ',sme_valid_min
      write(*,*)' sme       valid_max is ',sme_valid_max
      write(*,*)' rqf       FillValue is ',rqf_FillValue
      write(*,*)
   endif

   retrieval_flag_count   = 0
   moisture_fill_count    = 0
   moisture_min_count     = 0
   moisture_max_count     = 0
   temperature_fill_count = 0
   temperature_min_count  = 0
   temperature_max_count  = 0
   badlat_count           = 0
   badlon_count           = 0
   moisture_count         = 0
   temperature_count      = 0

   COUNTLOOP: do icount=1,counts

      ! Soil Moisture

      moisture_good = .true.

      if (moisture(icount) == moisture_FillValue(1))  then
          moisture_fill_count = moisture_fill_count + 1
          moisture_good = .false.
      endif
      if (moisture(icount) < moisture_valid_min(1)) then
          moisture_min_count = moisture_min_count + 1
          moisture_good = .false.
      endif
      if (moisture(icount) > moisture_valid_max(1)) then
          moisture_max_count = moisture_max_count + 1
          moisture_good = .false.
!        write(*,*)'observation ',icount, ' value = ',moisture(icount), ' > allowed max of ',moisture_valid_max(1)
      endif
      if (retrieval_flag(icount) == rqf_FillValue(1)) then
          retrieval_flag_count = retrieval_flag_count + 1
          moisture_good = .false.
      endif
!     if (moisture_error_std(icount) == sme_FillValue(1)) cycle COUNTLOOP
!     if (moisture_error_std(icount) < sme_valid_min(1)) cycle COUNTLOOP
!     if (moisture_error_std(icount) > sme_valid_max(1)) cycle COUNTLOOP

      ! Surface Temperature

      temperature_good = .false.  ! easiest way to exclude these given
                                  ! the NSIDC disclaimer above.

      if (temperature(icount) == temperature_FillValue(1))  then
          temperature_fill_count = temperature_fill_count + 1
          temperature_good = .false.
      endif
      if (temperature(icount) < temperature_valid_min(1)) then
          temperature_min_count = temperature_min_count + 1
          temperature_good = .false.
      endif
      if (temperature(icount) > temperature_valid_max(1)) then
          temperature_max_count = temperature_max_count + 1
          temperature_good = .false.
      endif

      ! Rejection section

      rlat = real(latitude( icount),r8)
      rlon = real(longitude(icount),r8)
      if (rlon < 0.0) rlon = rlon + 360.0_r8

      ! ensure the lat/longitude values are in range
      if ( rlat >  90.0_r8 .or. rlat <  -90.0_r8 ) then
         badlat_count = badlat_count + 1
         cycle COUNTLOOP
      endif
      if ( rlon > 360.0_r8 .or. rlon <    0.0_r8 ) then
         badlon_count = badlon_count + 1
         cycle COUNTLOOP
      endif

      call get_time(obs_time(icount), osec, oday)
      qc = retrieval_flag(icount)

      if ( moisture_good ) then

         ! Since the soil moisture error values are all _FillValue in the data
         ! I have, we are simply assuming something like 20% ...
         ! cannot have an error standard deviation of 0.0, so there must be an alternative
         ! minimum error specification ... just using 0.01 until proven otherwise
         ! Both 20% and 0.01 have no scientific basis and should be explored.

         obs_val = real(moisture(icount),r8)
         err_std = max(obs_val/20.0_r8, 0.01_r8)

         call create_3d_obs(rlat, rlon, sm_depth_m, VERTISHEIGHT, obs_val, &
                           SOIL_MOISTURE, err_std, oday, osec, qc, obs)

         call add_obs_to_seq(obs_seq, obs, obs_time(icount), prev_obs, prev_time, first_obs)

         moisture_count = moisture_count + 1
      endif

      if ( temperature_good ) then

         obs_val = real(temperature(icount),r8)
         err_std = max(obs_val*0.02_r8, 0.01_r8) ! pure speculation

         call create_3d_obs(rlat, rlon, t_depth_m, VERTISHEIGHT, obs_val, &
                           SOIL_TEMPERATURE, err_std, oday, osec, qc, obs)

         call add_obs_to_seq(obs_seq, obs, obs_time(icount), prev_obs, prev_time, first_obs)

         temperature_count = temperature_count + 1
      endif
   enddo COUNTLOOP

   if (verbose) then
      write(*,*)
      write(*,*)' retrieval_flag_count   is ', retrieval_flag_count
      write(*,*)' moisture_fill_count    is ', moisture_fill_count
      write(*,*)' moisture_min_count     is ', moisture_min_count
      write(*,*)' moisture_max_count     is ', moisture_max_count
      write(*,*)' temperature_fill_count is ', temperature_fill_count
      write(*,*)' temperature_min_count  is ', temperature_min_count
      write(*,*)' temperature_max_count  is ', temperature_max_count
      write(*,*)' badlat_count           is ', badlat_count
      write(*,*)' badlon_count           is ', badlon_count
      write(*,*)' leaves ', counts - moisture_fill_count - retrieval_flag_count - &
                    moisture_min_count - moisture_max_count - badlat_count - badlon_count
      write(*,*)' # GOOD soil_moisture   is ', moisture_count
      write(*,*)' # GOOD temperature     is ', temperature_count
      write(*,*)
   endif

   moisture_fill_count    = 0
   retrieval_flag_count   = 0
   moisture_min_count     = 0
   moisture_max_count     = 0
   temperature_fill_count = 0
   temperature_min_count  = 0
   temperature_max_count  = 0
   badlat_count           = 0
   badlon_count           = 0
   moisture_count         = 0
   temperature_count      = 0

   deallocate(dimlens)
   deallocate(longitude)
   deallocate(latitude)
   deallocate(temperature)
   deallocate(moisture)
   deallocate(moisture_error_std)
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


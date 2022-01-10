! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program LPRM_L3_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     LPRM_L3_to_obs - these data files were provided by Daniel Hagan, NUIST
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, MISSING_R8

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time, &
                              increment_time, get_time, set_date, operator(-),   &
                              print_date, operator(+)

use     utilities_mod, only : initialize_utilities, find_namelist_in_file,       &
                              check_namelist_read, nmlfileunit, do_output,       &
                              get_next_filename, error_handler, E_ERR, E_MSG,    &
                              find_textfile_dims, finalize_utilities,            &
                              do_nml_file, do_nml_term

use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, &
                                  nc_check, nc_get_dimension_size

use  obs_utilities_mod, only : getvar_real, get_2Dshort_as_r8, create_3d_obs

use      location_mod, only : VERTISHEIGHT, set_location

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,         &
                              static_init_obs_sequence, init_obs, destroy_obs,   &
                              write_obs_seq, init_obs_sequence, get_num_obs,     &
                              insert_obs_in_seq, destroy_obs_sequence,           &
                              set_copy_meta_data, set_qc_meta_data, set_qc,      & 
                              set_obs_values, set_obs_def, insert_obs_in_seq

use      obs_kind_mod, only : LPRM_SOIL_MOISTURE

use netcdf

implicit none

character(len=*), parameter :: source   = 'LPRM_L3_to_obs.f90'

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

character(len=256) :: next_infile
character(len=512) :: string1, string2

integer :: ncid, num_files, io, iunit, filenum
integer :: day, seconds
integer :: num_new_obs, obs_num, num_obs_in_1_file
integer :: i, j, k, num_Latitudes, num_Longitudes

logical :: first_obs, from_list = .false.

real(r8) :: oerr, qc, lat, lon, vertvalue
real(r8) :: smx_fill, smxe_fill

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time

real(r8), allocatable :: Latitude(:), Longitude(:)
real(r8), allocatable :: soil_moisture_x(:,:)
real(r8), allocatable :: sm_x_error(:,:)

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=256) :: input_file      = ''
character(len=256) :: input_file_list = ''
character(len=256) :: output_file     = 'obs_seq.lprm'
logical            :: debug           = .false.
real(r8)           :: lon_bounds(2)   = (/   0.0_r8, 360.0_r8 /)
real(r8)           :: lat_bounds(2)   = (/ -90.0_r8,  90.0_r8 /)

! representation error, which depends on model/data grid size, 
real(r8) :: lprm_rep_error = 0.1  ! representativeness error variance

namelist /LPRM_L3_to_obs_nml/ input_file, &
                          input_file_list, output_file, &
                          debug, lprm_rep_error, &
                          lon_bounds, lat_bounds

! start of executable code

call set_calendar_type(GREGORIAN)

!  read the necessary parameters from input.nml
call initialize_utilities()

call find_namelist_in_file('input.nml', 'LPRM_L3_to_obs_nml', iunit)
read(iunit, nml = LPRM_L3_to_obs_nml, iostat = io)
call check_namelist_read(iunit, io, 'LPRM_L3_to_obs_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=LPRM_L3_to_obs_nml)
if (do_nml_term()) write(     *     , nml=LPRM_L3_to_obs_nml)

! cannot have both a single filename and a list;
! the namelist must shut one off.

if (input_file /= '' .and. input_file_list /= '') then
  call error_handler(E_ERR, 'LPRM_L3_to_obs',                     &
                     'One of input_file or filelist must be NULL', &
                     source)
endif
if (input_file_list /= '') from_list = .true.

!-----------------------------------------------------------------------
! Get number of observations in a single file 

if (from_list) then
   next_infile = get_next_filename(input_file_list, 1)
   call find_textfile_dims(input_file_list, num_files)
else
   next_infile = input_file
   num_files   = 1
endif

ncid           = nc_open_file_readonly(next_infile, 'LPRM_L3_to_obs')
num_Latitudes  = nc_get_dimension_size(ncid,'Latitude')
num_Longitudes = nc_get_dimension_size(ncid,'Longitude')
call nc_close_file(ncid)

num_obs_in_1_file = num_Latitudes*num_Longitudes

num_new_obs = num_obs_in_1_file * num_files

call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)

do k = 1, num_copies
  call set_copy_meta_data(obs_seq, k, 'observation')
enddo

do k = 1, num_qc
  call set_qc_meta_data(obs_seq, k, 'original QC')
enddo

obs_num = 1
first_obs = .true.
   
! main loop that does either a single file or a list of files

filenum = 1
fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(input_file_list, filenum)
   else
      next_infile = input_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop

   ncid           = nc_open_file_readonly(next_infile, 'LPRM_L3_to_obs')
   num_Latitudes  = nc_get_dimension_size(ncid,'Latitude')
   num_Longitudes = nc_get_dimension_size(ncid,'Longitude')
 
   allocate(Latitude(       num_Latitudes))
   allocate(Longitude(      num_Longitudes))
   allocate(soil_moisture_x(num_Latitudes, num_Longitudes))
   allocate(sm_x_error(     num_Latitudes, num_Longitudes))

   call getvar_real(ncid, 'Longitude', Longitude)
   call getvar_real(ncid, 'Latitude',  Latitude)
   call get_2Dshort_as_r8(ncid, 'soil_moisture_x', soil_moisture_x, smx_fill)
   call get_2Dshort_as_r8(ncid, 'sm_x_error',      sm_x_error,      smxe_fill)
   call nc_close_file(ncid)

   ! Daniel wants to convert from 
   !   soil_moisture_x:long_name = "Volumetric Soil Moisture from X-band" ;
   !   soil_moisture_x:units = "percent" ;

   where( soil_moisture_x /= smx_fill ) soil_moisture_x = soil_moisture_x/100.0_r8

   ! time is determined from the file name
  
   obs_time = get_time_from_name(next_infile)
 
   call get_time(obs_time, seconds, day)
   
   ! ensure longitudes are [0,360] 
   where(Longitude < 0.0_r8) Longitude = Longitude + 360.0_r8

   LONLOOP: do i = 1, num_Longitudes
      lon = Longitude(i)
      if ( lon < lon_bounds(1) ) cycle LONLOOP
      if ( lon > lon_bounds(2) ) cycle LONLOOP

   LATLOOP: do j = 1, num_Latitudes
  
      lat = Latitude( j)

      if ( lat < lat_bounds(1) ) cycle LATLOOP
      if ( lat > lat_bounds(2) ) cycle LATLOOP

      qc  = 0.0

      if ( lat < lat_bounds(1) ) cycle LATLOOP

      !------------------------------------------------------------------
      ! observation error is instrument plus representativeness

      oerr  = sm_x_error(j,i) + lprm_rep_error
      vertvalue  = 0.005_r8  ! 5mm depth for X band

      if (soil_moisture_x(j,i) == smx_fill)  cycle LATLOOP
      if (     sm_x_error(j,i) == smxe_fill) cycle LATLOOP

      ! If you want to geographically subset, do it here.

      call create_3d_obs(lat, lon, vertvalue, VERTISHEIGHT, &
                soil_moisture_x(j,i), LPRM_SOIL_MOISTURE, &
                oerr, day, seconds, qc, obs)

      ! first one, insert with no prev.  otherwise, since all times are the
      ! same for this column, insert with the prev obs as the starting point.
      ! (the first insert with no prev means it will search for the right
      ! time ordered starting point.)
      if (first_obs) then
         call insert_obs_in_seq(obs_seq, obs)
         first_obs = .false.
      else
        call insert_obs_in_seq(obs_seq, obs, prev_obs)
      endif
      obs_num = obs_num+1
      prev_obs = obs

      if (obs_num > num_new_obs) then
         write(string1,*)'reached maximum number of observations while'
         write(string2,*)'reading file # ',filenum,' "'//trim(next_infile)//'"'
         call error_handler(E_ERR, 'LPRM_L3_to_obs', string1, &
                   source, text2=string2)
      endif

   enddo LATLOOP
   enddo LONLOOP

   filenum = filenum + 1

   deallocate(Latitude)
   deallocate(Longitude)
   deallocate(soil_moisture_x)
   deallocate(sm_x_error)

end do fileloop

! done with main loop.  if we added any obs to the sequence, write it out.
if (obs_num > 0) then
   print *, 'ready to write, nobs = ', get_num_obs(obs_seq)
   if (get_num_obs(obs_seq) > 0) &
      call write_obs_seq(obs_seq, output_file)

   call destroy_obs(obs)
   if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
else
   call error_handler(E_ERR,'LPRM_L3_to_obs','no obs converted', &
              source)
endif

call error_handler(E_MSG,'LPRM_L3_to_obs','Finished successfully.')
call finalize_utilities()


contains


!-----------------------------------------------------------------------
!> LPRM-TMI_L3_NT_SOILM3_V001-20120407T210538Z_20000724.nc.nc
!>                                             YYYYMMDD


function get_time_from_name(fname)
character(len=*), intent(in) :: fname
type(time_type) :: get_time_from_name

integer :: indx, indx1, indx2
integer :: year, month, day, hours, minutes, seconds

indx = index(fname,'.nc') ! may not be universally applicable

if (indx < 1) then
   call error_handler(E_ERR,'get_time_from_name','cannot determine where to start', &
              source, text2=fname)
endif

indx1 = indx - 8
indx2 = indx 

write(*,*)'subset of time string is ',fname(indx1:indx2)

read(fname(indx1:indx2),'(i4,2(i2))')year,month,day

! According to Daniel
hours   = 1
minutes = 30
seconds = 0

get_time_from_name = set_date(year,month,day,hours,minutes,seconds)

end function get_time_from_name


end program LPRM_L3_to_obs


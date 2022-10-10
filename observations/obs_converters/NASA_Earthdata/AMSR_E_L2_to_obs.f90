! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program AMSR_E_L2_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     AMSR_E_L2_to_obs - 
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

use  obs_utilities_mod, only : get_2Dshort_as_r8, create_3d_obs

use      location_mod, only : VERTISHEIGHT, set_location

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,         &
                              static_init_obs_sequence, init_obs, destroy_obs,   &
                              write_obs_seq, init_obs_sequence, get_num_obs,     &
                              insert_obs_in_seq, destroy_obs_sequence,           &
                              set_copy_meta_data, set_qc_meta_data, set_qc,      & 
                              set_obs_values, set_obs_def, insert_obs_in_seq

use      obs_kind_mod, only : AMSRE_A_SOIL_MOISTURE_X, AMSRE_D_SOIL_MOISTURE_X, &
                              AMSRE_A_SOIL_MOISTURE_C, AMSRE_D_SOIL_MOISTURE_C

use netcdf

implicit none

character(len=*), parameter :: source   = 'AMSR_E_L2_to_obs.f90'

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

character(len=256) :: next_infile
character(len=512) :: string1, string2

integer :: ncid, num_files, io, iunit, filenum
integer :: day, seconds
integer :: num_new_obs, obs_num, num_obs_in_1_file
integer :: i, j, k, nscan, npix
integer :: c_observation_kind, x_observation_kind

logical :: first_obs, from_list = .false.

real(r8) :: oerr, qc, lat, lon, vertvalue
real(r8) :: odc_fill, odx_fill, smc_fill, smce_fill, smx_fill, smxe_fill

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time

real(r8), allocatable :: Latitude(:,:), Longitude(:,:)
real(r8), allocatable :: rfi_code(:,:)
real(r8), allocatable :: soil_moisture_c(:,:)
real(r8), allocatable :: soil_moisture_x(:,:)
real(r8), allocatable :: sm_c_error(:,:)
real(r8), allocatable :: sm_x_error(:,:)
real(r8), allocatable :: opt_depth_c(:,:)
real(r8), allocatable :: opt_depth_x(:,:)

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

character(len=256) :: input_file      = ''
character(len=256) :: input_file_list = ''
character(len=256) :: output_file     = 'obs_seq.amsre'
logical            :: debug           = .false.
integer            :: max_rfi_code    = 2

! representation error, which depends on model/data grid size, 
real(r8) :: amsre_rep_error = 0.3  ! representativeness error, in percent

namelist /AMSR_E_L2_to_obs_nml/ input_file, &
                          input_file_list, output_file, &
                          debug, amsre_rep_error, max_rfi_code

! start of executable code

call set_calendar_type(GREGORIAN)

!  read the necessary parameters from input.nml
call initialize_utilities()

call find_namelist_in_file('input.nml', 'AMSR_E_L2_to_obs_nml', iunit)
read(iunit, nml = AMSR_E_L2_to_obs_nml, iostat = io)
call check_namelist_read(iunit, io, 'AMSR_E_L2_to_obs_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=AMSR_E_L2_to_obs_nml)
if (do_nml_term()) write(     *     , nml=AMSR_E_L2_to_obs_nml)

! cannot have both a single filename and a list;
! the namelist must shut one off.

if (input_file /= '' .and. input_file_list /= '') then
  call error_handler(E_ERR, 'AMSR_E_L2_to_obs',                     &
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

ncid  = nc_open_file_readonly(next_infile, 'AMSR_E_L2_to_obs')
nscan = nc_get_dimension_size(ncid,'nscan')
npix  = nc_get_dimension_size(ncid,'npix')
call nc_close_file(ncid)

num_obs_in_1_file = 2*nscan*npix*1.1 ! a little extra

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

   ! Need to know if it is ascending or descending 
   c_observation_kind = asc_desc(next_infile,'C') ! AMSRE_[A,D]_SOIL_MOISTURE_C
   x_observation_kind = asc_desc(next_infile,'X') ! AMSRE_[A,D]_SOIL_MOISTURE_X

   ncid  = nc_open_file_readonly(next_infile, 'AMSR_E_L2_to_obs')
   nscan = nc_get_dimension_size(ncid,'nscan')
   npix  = nc_get_dimension_size(ncid,'npix')
 
   allocate(Latitude(       npix,nscan))
   allocate(Longitude(      npix,nscan))
   allocate(rfi_code(       npix,nscan))
   allocate(soil_moisture_c(npix,nscan))
   allocate(soil_moisture_x(npix,nscan))
   allocate(sm_c_error(     npix,nscan))
   allocate(sm_x_error(     npix,nscan))
   allocate(opt_depth_c(    npix,nscan))
   allocate(opt_depth_x(    npix,nscan))

   call get_2Dshort_as_r8(ncid, 'Longitude', Longitude)
   call get_2Dshort_as_r8(ncid, 'Latitude',  Latitude)
   call get_2Dshort_as_r8(ncid, 'rfi_code',  rfi_code)
   call get_2Dshort_as_r8(ncid, 'soil_moisture_c', soil_moisture_c, smc_fill)
   call get_2Dshort_as_r8(ncid, 'soil_moisture_x', soil_moisture_x, smx_fill)
   call get_2Dshort_as_r8(ncid, 'sm_c_error',      sm_c_error,      smce_fill)
   call get_2Dshort_as_r8(ncid, 'sm_x_error',      sm_x_error,      smxe_fill)
   call get_2Dshort_as_r8(ncid, 'opt_depth_c',     opt_depth_c,     odc_fill)
   call get_2Dshort_as_r8(ncid, 'opt_depth_x',     opt_depth_x,     odx_fill)
   call nc_close_file(ncid)

   ! time is determined from the file name
  
   obs_time = get_time_from_name(next_infile)
 
   call get_time(obs_time, seconds, day)
   
   ! ensure longitudes are [0,360] 
   where(Longitude < 0.0_r8) Longitude = Longitude + 360.0_r8

   scanloop: do j = 1, nscan
   pixloop: do i = 1, npix
  
      if (rfi_code(i,j) > max_rfi_code ) cycle pixloop

      lon = Longitude(i,j)
      lat = Latitude( i,j)
      qc  = rfi_code( i,j)

      !------------------------------------------------------------------
      ! Create C-band observation
      ! observation error is instrument plus representativeness

      oerr  = sm_c_error(i,j) + amsre_rep_error  ! units of standard deviation
      vertvalue  = 0.005_r8  ! 5mm depth for C band

      if ((soil_moisture_c(i,j) /= smc_fill) .and. &
               (sm_c_error(i,j) /= smce_fill) ) then
           call create_3d_obs(lat, lon, vertvalue, VERTISHEIGHT, &
                   soil_moisture_c(i,j), c_observation_kind, &
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
            call error_handler(E_ERR, 'AMSR_E_L2_to_obs', string1, &
                       source, text2=string2)
         endif

      endif

      !------------------------------------------------------------------
      ! Create X-band observation
      ! observation error is instrument plus representativeness

      oerr  = sm_x_error(i,j) + amsre_rep_error
      vertvalue  = 0.005_r8  ! 5mm depth for X band

      if ((soil_moisture_x(i,j) /= smx_fill) .and. &
               (sm_x_error(i,j) /= smxe_fill) ) then
           call create_3d_obs(lat, lon, vertvalue, VERTISHEIGHT, &
                   soil_moisture_x(i,j), x_observation_kind, &
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
            call error_handler(E_ERR, 'AMSR_E_L2_to_obs', string1, &
                      source, text2=string2)
         endif

      endif

   enddo pixloop
   enddo scanloop

   filenum = filenum + 1

   deallocate(Latitude)
   deallocate(Longitude)
   deallocate(rfi_code)
   deallocate(soil_moisture_c)
   deallocate(soil_moisture_x)
   deallocate(sm_c_error)
   deallocate(sm_x_error)
   deallocate(opt_depth_c)
   deallocate(opt_depth_x)

end do fileloop

! done with main loop.  if we added any obs to the sequence, write it out.
if (obs_num > 0) then
   print *, 'ready to write, nobs = ', get_num_obs(obs_seq)
   if (get_num_obs(obs_seq) > 0) &
      call write_obs_seq(obs_seq, output_file)

   call destroy_obs(obs)
   if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
else
   call error_handler(E_ERR,'AMSR_E_L2_to_obs','no obs converted', &
              source)
endif

call error_handler(E_MSG,'AMSR_E_L2_to_obs','Finished successfully.')
call finalize_utilities()


contains 

!-----------------------------------------------------------------------
!> determine the observation kind from the filename. There are separate
!> kinds for ascending vs. descending orbits.
!>
!> LPRM-AMSR_E_L2_A_SOILM2_V002_20030702114630.nc
!> LPRM-AMSR_E_L2_D_SOILM2_V002_20030702114630.nc

function asc_desc(fname,band)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: band
integer :: asc_desc

character(len=64) :: mystring

if (index(fname,'L2_A_') > 0 ) then
   mystring = 'AMSRE_A_SOIL_MOISTURE_'//trim(band)
elseif (index(fname,'L2_D_') > 0 ) then
   mystring = 'AMSRE_D_SOIL_MOISTURE_'//trim(band)
else
   call error_handler(E_ERR,'asc_desc','cannot find asc/desc character', &
              source, text2=fname)
endif

SELECT CASE (mystring)
   CASE        ('AMSRE_A_SOIL_MOISTURE_C')
      asc_desc = AMSRE_A_SOIL_MOISTURE_C
   CASE        ('AMSRE_A_SOIL_MOISTURE_X')
      asc_desc = AMSRE_A_SOIL_MOISTURE_X
   CASE        ('AMSRE_D_SOIL_MOISTURE_C')
      asc_desc = AMSRE_D_SOIL_MOISTURE_C
   CASE        ('AMSRE_D_SOIL_MOISTURE_X')
      asc_desc = AMSRE_D_SOIL_MOISTURE_X
   CASE DEFAULT
      call error_handler(E_ERR,'asc_desc','bad string "'//trim(mystring)//'', &
                 source)
END SELECT

end function asc_desc


!-----------------------------------------------------------------------
!> determine the observation kind from the filename. There are separate
!> kinds for ascending vs. descending orbits.
!>
!> LPRM-AMSR_E_L2_D_SOILM2_V002_20030702114630.nc
!>                              YYYYMMDDHHmmss

function get_time_from_name(fname)
character(len=*), intent(in) :: fname
type(time_type) :: get_time_from_name

integer :: indx, indx1, indx2
integer :: year, month, day, hours, minutes, seconds

indx = index(fname,'_V002_')

if (indx < 1) then
   call error_handler(E_ERR,'get_time_from_name','cannot determine where to start', &
              source, text2=fname)
endif

indx1 = indx + 6 ! need to know where '_V002_' ends
indx2 = indx1 + 13 

read(fname(indx1:indx2),'(i4,5(i2))')year,month,day,hours,minutes,seconds

get_time_from_name = set_date(year,month,day,hours,minutes,seconds)

end function get_time_from_name



end program AMSR_E_L2_to_obs


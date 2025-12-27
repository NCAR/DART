! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! This converter supports in-situ ocean observations for the InaCAWO project.
! The types of data that can be ingested are:  
! 1. ARVOR-C: Argo-type profiling floats (T/S profiles)
! 2. ARVOR-I: Argo-type (ice-capable) floats (T/S profiles)
!
! The raw data files come in ASCII format, containing 
! a header with float metadata and multiple rows representing
! profile levels during a dive/ascent cycle
! 
! The raw data comes in with several characteristics such as timestamp, 
! instrument ID information, geographic location, cycle metadata (e.g., 
! start time, mission settings, descent slices, parking depth, etc).
! The converter typically extracts: time, lat, lon, pressure/depth, 
! temperature, and salinity. 
!
! Note on units:  
! temp in Kelvin: convert to C.
! pres in dbar: convert to depth (m, +down).

program arvor_to_obs

use types_mod,            only : r8, t_kelvin, MISSING_R8
use time_manager_mod,     only : time_type, set_calendar_type, GREGORIAN, get_time,   &   
                                 set_date, print_date, operator(+), operator(-)
use utilities_mod,        only : initialize_utilities, find_namelist_in_file,         &   
                                 nmlfileunit, error_handler, do_nml_term, E_ERR,      &
                                 finalize_utilities, do_nml_file, get_next_filename,  &
                                 find_textfile_dims, file_exist, E_MSG
use location_mod,         only : VERTISHEIGHT, set_location
use obs_sequence_mod,     only : obs_type, obs_sequence_type, init_obs, get_num_obs,  &
                                 static_init_obs_sequence, init_obs_sequence,         &
                                 set_copy_meta_data, set_qc_meta_data, write_obs_seq, &
                                 destroy_obs_sequence, insert_obs_in_seq,             &
                                 set_obs_values, set_obs_def, destroy_obs
use obs_utilities_mod,    only : create_3d_obs, add_obs_to_seq
use obs_kind_mod,         only : FLOAT_TEMPERATURE, FLOAT_SALINITY
use obs_def_mod,          only : obs_def_type, set_obs_def_time, set_obs_def_key,     &
                                 set_obs_def_error_variance, set_obs_def_location,    &
                                 set_obs_def_type_of_obs
use read_csv_mod,         only : csv_file_type, csv_get_nrows, csv_get_field,         &
                                 csv_open, csv_close, csv_print_header

implicit none

character(len=*), parameter :: source = 'arvor_to_obs'

character(len=512) :: string1, string2

! File variables
character(len=256) :: next_infile
integer            :: io, iunit, filenum
integer            :: num_new_obs, nfiles
integer            :: num_valid_obs
logical            :: from_list = .false., first_obs = .true.

! Obs sequence
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time
integer                 :: num_copies = 1, &     ! number of copies in sequence
                           num_qc     = 1        ! number of QC entries
real(r8)                :: obs_qc

! Data arrays
real(r8), allocatable   :: lon(:), lat(:), pres(:) ! Location 
real(r8), allocatable   :: temp(:), saln(:)        ! Physical observations

character(len=512), allocatable :: dat(:)

integer,  dimension(2)  :: float_types = [FLOAT_TEMPERATURE, FLOAT_SALINITY]
real(r8), dimension(2)  :: float_errs  

! ARVOR profile
type arvor
    type(time_type) :: dat 
    real(r8) :: lat
    real(r8) :: lon
    real(r8) :: dep
    real(r8) :: tmp
    real(r8) :: sal 
end type arvor

type(arvor), allocatable :: prof(:)

! csv obs file 
type(csv_file_type) :: cf

!------------------------------------------------------------------------
!  Declare namelist parameters
character(len=256) :: file_in           = ''
character(len=256) :: file_list         = ''
character(len=256) :: file_out          = 'obs_seq.arvor'
integer            :: avg_obs_per_file  = 500000
real(r8)           :: obs_error_temp    = 0.02_r8  ! Profiling floats are highly accurate
real(r8)           :: obs_error_sal     = 0.02_r8  ! ARVOR CTDs are precise
logical            :: debug             = .true.


namelist /arvor_to_obs_nml/ file_in,          &
                            file_list,        &
                            obs_error_temp,   &
                            obs_error_sal,    &
                            file_out,         &
                            avg_obs_per_file, &
                            debug

! Start Converter
call initialize_utilities()

! Read the namelist options
call find_namelist_in_file('input.nml', 'arvor_to_obs_nml', iunit)
read(iunit, nml = arvor_to_obs_nml, iostat = io)

if (do_nml_file()) write(nmlfileunit, nml=arvor_to_obs_nml)
if (do_nml_term()) write(     *     , nml=arvor_to_obs_nml)

! Set the calendar kind
call set_calendar_type(GREGORIAN)

! Check the files
if (file_in   /= '' .and. file_list /= '') then
   string1 = 'One of input ARVOR file or the file list must be NULL'
   call error_handler(E_ERR, source, string1)
endif
if (file_list /= '') from_list = .true.

! Get number of observations
num_new_obs = avg_obs_per_file
if (from_list) then
   call find_textfile_dims(file_list, nfiles)
   num_new_obs = avg_obs_per_file * nfiles
endif

! Assign obs errors
float_errs = [obs_error_temp, obs_error_sal]

! Initialize obs seq file
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

print *
if (file_exist(file_out)) then
   write(*, '(A)') 'Output file: '//trim(adjustl(file_out))//' exists. Replacing it ...'
else
   write(*, '(A)') 'Creating "'//trim(adjustl(file_out))//'" file.'
endif

call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
call set_copy_meta_data(obs_seq, num_copies, 'ARVOR observation')
call set_qc_meta_data(obs_seq, num_qc, 'ARVOR QC')

! Loop over the obs files
filenum = 1
obs_qc  = 0.0_r8

FILELOOP: do
  ! Get the data file
  if (from_list) then
     next_infile = get_next_filename(file_list, filenum)
  else
     next_infile = file_in
     if (filenum > 1) next_infile = ''
  endif
  if (next_infile == '') exit FILELOOP
  print *
  if (debug) write(*, '(A, i0, X, A)') 'Input file: #', filenum, trim(next_infile)

  ! Read data from file
  call get_profile_data(next_infile)  

  ! Add obs to sequence
  call fill_obs(obs, prev_obs, prev_time, obs_qc)

  call cleanup

  filenum = filenum + 1
enddo FILELOOP

call finish_obs()
call finalize_utilities()


contains


!-----------------------------------------------
! Read the data profile and the associated dates 
subroutine get_profile_data(filename)

character(len=*), intent(in) :: filename
character(len=*), parameter  :: routine = 'get_profile_data'

integer :: i, k, nobs

! Open csv file and get dims
call csv_open(filename, cf, context=routine)
nobs = csv_get_nrows(cf)

if (debug) call csv_print_header(cf)

allocate(dat(nobs) , prof(nobs))
allocate(lat(nobs) , lon(nobs), pres(nobs))
allocate(temp(nobs), saln(nobs))

! Read the data
call csv_get_field(cf, 'lat', lat, routine)
call csv_get_field(cf, 'lon', lon, routine)
call csv_get_field(cf, 'sea_water_pressure', pres, routine)
call csv_get_field(cf, 'sea_water_temperature', temp, routine)
call csv_get_field(cf, 'salinity', saln, routine)
call csv_get_field(cf, 'dat', dat, routine)

call csv_close(cf)

! Data adjustments
k = 0
do i = 1, nobs
   if (lat(i)  == MISSING_R8 .or. &
       lon(i)  == MISSING_R8 .or. &
       pres(i) == MISSING_R8 .or. &
       temp(i) == MISSING_R8 .or. & 
       saln(i) == MISSING_R8 ) cycle
   k = k + 1

   if (lon(i) < 0.0_r8) lon(i) = lon(i) + 360.0_r8

   prof(k)%lon = lon(i)
   prof(k)%lat = lat(i)
   prof(k)%dat = parse_time(dat(i))
   prof(k)%dep = depth_from_pressure(pres(i), lat(i))
   prof(k)%tmp = temp(i) - t_kelvin
   prof(k)%sal = saln(i) ! Salinity is in PSU; consistent with ROMS
enddo
num_valid_obs = k

if (debug) then
   write(*, '(A, i0)') '  >> Total number of obs: ', nobs
   write(*, '(A, i0)') '  >> Valid observations count: ', num_valid_obs  

   do i = 1, num_valid_obs
      write(string1, '(5X, A, i3, 5(A,f10.4), A)') '* obs #', i, ', lat:', prof(i)%lat, &
                                          ', lon:', prof(i)%lon, ', dep:', prof(i)%dep, &
                                          ', T:'  , prof(i)%tmp, ', S:'  , prof(i)%sal, &
                                          ', date: '
      call print_date(prof(i)%dat, str=string1)
   enddo
endif 

end subroutine get_profile_data


!------------------------------------------------------
! Given a profile, add the metadata in the obs sequence  
subroutine fill_obs(obs, prev_obs, prev_time, oqc)

type(obs_type),        intent(inout) :: obs, prev_obs
type(time_type),       intent(inout) :: prev_time
real(r8),              intent(in)    :: oqc

type(time_type) :: odat
real(r8)        :: olon, olat, odep, oval, oerr
integer         :: otype, iobs, osec, oday, k

! T, S
do k = 1, 2
   otype = float_types(k)
   oerr  = float_errs(k)

   do iobs = 1, num_valid_obs
      olon = prof(iobs)%lon
      olat = prof(iobs)%lat
      odep = prof(iobs)%dep
      odat = prof(iobs)%dat

      if (k == 1) oval = prof(iobs)%tmp
      if (k == 2) oval = prof(iobs)%sal

      call get_time(odat, osec, oday)
      call create_3d_obs(olat, olon, odep, VERTISHEIGHT, oval, otype, oerr, oday, osec, oqc, obs)
      call add_obs_to_seq(obs_seq, obs, odat, prev_obs, prev_time, first_obs)
   enddo
enddo

end subroutine fill_obs


!------------------------------------------------------------
! All collected obs (if any) are written in the seq file. 
subroutine finish_obs()

integer :: obs_num

obs_num = get_num_obs(obs_seq)

if (obs_num > 0) then
   print *
   if (debug) write(*, '(A, i0, A)') '> Ready to write ', obs_num, ' observations:'

   call write_obs_seq(obs_seq, file_out)
   call destroy_obs(obs)
   call destroy_obs_sequence(obs_seq)
else
   string1 = 'No obs were converted.'
   call error_handler(E_ERR, source, string1)
endif
call error_handler(E_MSG, source, 'Finished successfully.')

end subroutine finish_obs


!-----------------------------------------
! Using the date string, set the DART time
function parse_time(datestr) result(obs_time)

character(len=*), parameter :: routine = 'parse_time'

character(len=*), intent(in) :: datestr 
type(time_type)              :: obs_time

integer :: year, month, day, hour, minute, second

!example: 2025-10-06T06:10:00Z

read(datestr, '(i4, 5(1x,i2))', iostat=io) year, month, day, hour, minute, second
if (io /= 0) then
   write(string1, *) 'Unable to read time variable. Error status was ', io
   call error_handler(E_ERR, routine, string1, source)
endif

obs_time = set_date(year, month, day, hour, minute, second)

end function parse_time


!--------------------------------------------------------------
! Convert sea water pressure to depth (+ve, downward) in meters
! p_dbar : pressure in decibar (dbar)
! lat_deg: latitude in degrees ([-90, 90])
!
! Includes latitude-dependent gravity and the 1.092e-6*p compressibility term.
! Valid for typical oceanic ranges (0–12,000 dbar).
function depth_from_pressure(p_dbar, lat_deg) result(depth_m)

character(len=*), parameter :: routine = 'depth_from_pressure'

real(r8), intent(in) :: p_dbar
real(r8), intent(in) :: lat_deg
real(r8)             :: depth_m
real(r8)             :: x, sin2, sin4, g

if (p_dbar < 0._r8 .or. abs(lat_deg) > 90._r8) then
   write(string1, *) 'Cannot compute depth with input pressure: ', p_dbar 
   write(string2, *) 'and latitude: ', lat_deg 
   call error_handler(E_ERR, routine, string1, source, text2 = string2)
end if

! Latitude-dependent gravity (m/s^2), Saunders (1981)/UNESCO
x    = sin(lat_deg * (acos(-1._r8)/180._r8))     ! sin(phi)
sin2 = x*x
sin4 = sin2*sin2
g    = 9.780318_r8 * (1._r8 + 5.2788e-3_r8*sin2 + 2.36e-5_r8*sin4) + &
       1.092e-6_r8 * p_dbar

! Polynomial (UNESCO 1983): depth (m) from pressure (dbar)
depth_m = (((-1.82e-15_r8*p_dbar + 2.279e-10_r8)*p_dbar - 2.2512e-5_r8)*p_dbar + 9.72659_r8) &
          * p_dbar / g

end function depth_from_pressure


!------------------------------------------------------------
! Clear up memory 
subroutine cleanup()

if (allocated(prof)) deallocate(prof)
if (allocated(lon) ) deallocate(lon)
if (allocated(lat) ) deallocate(lat)
if (allocated(dat) ) deallocate(dat)
if (allocated(temp)) deallocate(temp)
if (allocated(saln)) deallocate(saln)
if (allocated(pres)) deallocate(pres)

end subroutine cleanup

end program arvor_to_obs

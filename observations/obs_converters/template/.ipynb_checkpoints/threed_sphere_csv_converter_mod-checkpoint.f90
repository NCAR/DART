! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program threed_csv_converter_mod

!> title = "Generalized 3D Sphere CSV Observation Converter"
!> institution = "NCAR" ;
!> source = "NCAR/DAReS" ;
!> comment = "Generalized converter for 3D sphere data from CSV files" ;
!> references = "http://www.image.ucar.edu/DAReS/DART/DART_download" ;
!> dataset_title = "Generalized 3D Sphere CSV Data" ;

use         types_mod, only : r8, digits12, MISSING_R8

use     utilities_mod, only : initialize_utilities, find_namelist_in_file, &
                              check_namelist_read, nmlfileunit, &
                              error_handler, E_ERR, E_MSG, &
                              finalize_utilities, do_nml_file, do_nml_term, &
                              open_file, close_file

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_time, get_time, print_time, &
                              set_date, get_date, print_date, &
                              operator(+), operator(-)
                              
use      location_mod, only : get_location, location_type, set_location, VERTISSURFACE, VERTISHEIGHT, VERTISUNDEF

use  obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs, &
                              static_init_obs_sequence, init_obs_sequence, &
                              set_copy_meta_data, set_qc_meta_data, &
                              get_num_obs, write_obs_seq, destroy_obs_sequence
                              
use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs

use obs_kind_mod, only : YOUR_OBS_KIND, get_kind_index

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = & '$URL$'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'
character(len=*), parameter :: routine  = 'csv_to_obs'

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

real(r8), parameter :: qc = 0.0_r8        ! default QC value

! variables for file handling and observations
character(len=256) :: output_file
character(len=512) :: string1, string2

integer :: iunit, io, nlines, iline
integer :: oday, osec, year, month, day, hour, minutes, seconds
integer :: num_new_obs, nmissing, subsample_intv, ncolumns
integer :: lat_col, lon_col, time_col, var_col, qc_col

logical :: first_obs

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, prev_time, base_time, delta_time

real(digits12), allocatable :: time(:)
real(r8), allocatable :: lat(:), lon(:)
real(r8), allocatable :: your_variable(:)
real(r8) :: missing_value, your_variable_error_std

!------------------------------------------------------------------------
!  Declare namelist parameters

character(len=256) :: input_file       = 'input.csv'
character(len=256) :: output_file_base = 'obs_seq'
logical            :: debug            = .false.

namelist /csv_to_obs_nml/ input_file, output_file_base, &
                          debug, subsample_intv, your_variable_error_std

!------------------------------------------------------------------------
! Start of executable code
!------------------------------------------------------------------------

! read and record standard parameters from input namelist
call initialize_utilities('csv_to_obs', .true., .true.)
call find_namelist_in_file('input.nml', 'csv_to_obs_nml', iunit)
read(iunit, nml = csv_to_obs_nml, iostat = io)

if (do_nml_file()) write(nmlfileunit, nml=csv_to_obs_nml)
if (do_nml_term()) write(     *     , nml=csv_to_obs_nml)

call set_calendar_type(GREGORIAN)

! open the CSV file and read number of lines
iunit = open_file(input_file, 'formatted', 'read')
nlines = count_file_lines(iunit)

! allocate arrays for data
allocate(time(nlines), lat(nlines), lon(nlines), your_variable(nlines))

! initialize observation sequence
call static_init_obs_sequence()

! read the CSV file
call read_csv_file(iunit, nlines, time, lat, lon, your_variable, ncolumns, lat_col, lon_col, time_col, var_col, qc_col)

! close the CSV file
call close_file(iunit)

num_new_obs = nlines

! initialize observation sequence
call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
call set_copy_meta_data(obs_seq, 1, 'YOUR_VARIABLE observation')
call set_qc_meta_data(obs_seq, 1, 'YOUR_VARIABLE QC')

first_obs = .true.
nmissing = 0
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! loop over each row in the CSV file
obsloop: do iline = 1, nlines

   if (your_variable(iline) == missing_value) then
      nmissing = nmissing + 1
      cycle obsloop
   endif

   ! set observation location
   call set_location(location, lat(iline), lon(iline), VERTISUNDEF)
   ! create and add observation to sequence
   call create_3d_obs(location, your_variable(iline), VERTISUNDEF, obs_time, obs_seq, first_obs)

end do obsloop

! if we added obs to the sequence, write it out to a file
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   if (debug) print *, '                  skipping = ', nmissing
   call write_obs_seq(obs_seq, output_file)
else
   write(string1,*)'no observations for output file'
   write(string2,*)'"'//trim(output_file)//'"'
   call error_handler(E_MSG, routine, string1, text2=string2)
endif

! destroy current observation sequence to prepare for the next time step
call destroy_obs_sequence(obs_seq)

call error_handler(E_MSG, routine, 'Finished successfully.')
call finalize_utilities()

contains

subroutine read_csv_file(iunit, nlines, time, lat, lon, var, ncolumns, lat_col, lon_col, time_col, var_col, qc_col)
  implicit none
  integer, intent(in) :: iunit, nlines, ncolumns
  integer, intent(out) :: lat_col, lon_col, time_col, var_col, qc_col
  real(digits12), intent(out) :: time(nlines), lat(nlines), lon(nlines)
  real(r8), intent(out) :: var(nlines)
  integer :: iline, ios
  character(len=256) :: line

  ! Determine the column indices from the header
  read(iunit, '(A)', iostat=ios) line
  if (ios /= 0) call error_handler(E_ERR, 'read_csv_file', 'Error reading CSV header')

  call parse_header(line, ncolumns, lat_col, lon_col, time_col, var_col, qc_col)

  ! Read the data
  do iline = 1, nlines
    read(iunit, *, iostat=ios) time(iline), lat(iline), lon(iline), var(iline)
    if (ios /= 0) then
      write(string1,*) 'Error reading CSV data at line ', iline
      call error_handler(E_ERR, 'read_csv_file', string1)
    endif
  end do
end subroutine read_csv_file

subroutine parse_header(line, ncolumns, lat_col, lon_col, time_col, var_col, qc_col)
  implicit none
  character(len=256), intent(in) :: line
  integer, intent(out) :: ncolumns, lat_col, lon_col, time_col, var_col, qc_col
  integer :: i, col_pos
  character(len=32) :: col_name
  character(len=1) :: delimiter
  logical :: found

  delimiter = ','

  lat_col = -1
  lon_col = -1
  time_col = -1
  var_col = -1
  qc_col = -1

  col_pos = 1
  ncolumns = 0

  do i = 1, len_trim(line)
    if (line(i:i) == delimiter .or. i == len_trim(line)) then
      if (i == len_trim(line)) col_name = trim(line(col_pos:i))
      if (trim(col_name) == 'lat' .or. trim(col_name) == 'latitude') lat_col = ncolumns + 1
      if (trim(col_name) == 'lon' .or. trim(col_name) == 'longitude') lon_col = ncolumns + 1
      if (trim(col_name) == 'time') time_col = ncolumns + 1
      if (trim(col_name) == 'var' .or. trim(col_name) == 'variable') var_col = ncolumns + 1
      if (trim(col_name) == 'qc') qc_col = ncolumns + 1

      col_pos = i + 1
      col_name = ''
      ncolumns = ncolumns + 1
    else
      col_name = trim(col_name) // line(i:i)
    endif
  end do

  if (lat_col == -1 .or. lon_col == -1 .or. time_col == -1 .or. var_col == -1 .or. qc_col == -1) then
    call error_handler(E_ERR, 'parse_header', 'Required columns not found in CSV header')
  endif
end subroutine parse_header

function count_file_lines(iunit) result(count)
  implicit none
  integer, intent(in) :: iunit
  integer :: count
  integer :: ios
  character(len=256) :: line

  count = 0
  do
    read(iunit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    count = count + 1
  end do

  rewind(iunit)
end function count_file_lines

subroutine create_3d_obs(location, value, vert_coord, obs_time, obs_seq, first_obs)

  use obs_def_mod, only : obs_def_type, set_obs_def_location, set_obs_def_type_of_obs, &
                          set_obs_def_time, set_obs_def_error_variance, set_obs_def
  use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
  use time_manager_mod, only : time_type, set_time
  use location_mod, only : location_type, set_location

  type(location_type), intent(in) :: location
  real(r8), intent(in) :: value
  integer, intent(in) :: vert_coord
  type(time_type), intent(in) :: obs_time
  type(obs_sequence_type), intent(inout) :: obs_seq
  logical, intent(inout) :: first_obs

  type(obs_def_type) :: obs_def
  type(obs_type) :: obs
  integer :: kind

  ! get the kind index for the variable of interest
  call get_kind_index('YOUR_VARIABLE_KIND', kind)
  call obs_def%initialize()
  call obs_def%set(location, kind, value)

  ! set observation time
  call obs_def%set_time(obs_time)

  ! add observation to sequence
  call add_obs_to_seq(obs_seq, obs, obs_time, first_obs)
  first_obs = .false.

end subroutine create_3d_obs

subroutine add_obs_to_seq(seq, obs, obs_time, first_obs)
  use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq
  use time_manager_mod, only : time_type, operator(>=)

  type(obs_sequence_type), intent(inout) :: seq
  type(obs_type),          intent(inout) :: obs
  type(time_type),         intent(in)    :: obs_time
  logical,                 intent(inout) :: first_obs

  if(first_obs) then    ! for the first observation, no prev_obs
    call insert_obs_in_seq(seq, obs)
    first_obs = .false.
  else
    call insert_obs_in_seq(seq, obs)
  endif

end subroutine add_obs_to_seq

end program threed_csv_converter_mod

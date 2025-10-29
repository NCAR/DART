! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! This converter supports in-situ ocean observations for the InaCAWO project.
! The types of data that can be ingested are: 
! 
! 1. SVP    : Surface Velocity Program difters (surface currents + SST)
! 2. ARVOR-C: Argo-type profiling floats (T/S profiles)
! 3. ARVOR-I: Argo-type (ice-capable) floats (T/S profiles)
!
! The raw data files for the 3 products come in ASCII format, containing 
! a header with float/drifter metadata and multiple rows representing: 
! - profile levels during a dive/ascent cycle, or
! - surface ocean conditions at different times. 
! 
! The raw data comes in with several characteristics such as timestamp, 
! instrument ID information, geographic location, cycle metadata (e.g., 
! start time, mission settings, descent slices, parking depth, etc).
! The converter typically extracts: time, lat, lon, pressure/depth, 
! temperature, salinity, U, and V. 
!
! Note on units:  
! ARVOR temp in Kelvin: convert to °C.
! ARVOR pres in dbar: convert to depth (m, +down).
!  SVP  sst in Kelvin: °C; u/v in m/s


program bmkg_to_obs

use types_mod,            only : r8, i8, i4, digits12, t_kelvin, MISSING_R8
use time_manager_mod,     only : time_type, set_calendar_type, GREGORIAN, set_time,      &   
                                 get_time, print_time, set_date, get_date, print_date,   &   
                                 write_time, operator(+), operator(-)
use utilities_mod,        only : initialize_utilities, find_namelist_in_file, E_ERR,     &   
                                 check_namelist_read, nmlfileunit, error_handler, E_MSG, &
                                 finalize_utilities, do_nml_file, get_next_filename,     &   
                                 do_nml_term, find_textfile_dims, file_exist, open_file, &
                                 close_file, to_upper, string_to_real
use location_mod,         only : VERTISSURFACE, VERTISHEIGHT, set_location
use obs_sequence_mod,     only : obs_type, obs_sequence_type, init_obs, get_num_obs,     &
                                 static_init_obs_sequence, init_obs_sequence,            &
                                 set_copy_meta_data, set_qc_meta_data, write_obs_seq,    &
                                 destroy_obs_sequence, insert_obs_in_seq, set_qc,        &
                                 set_obs_values, set_obs_def, destroy_obs
use obs_utilities_mod,    only : create_3d_obs, add_obs_to_seq
use obs_kind_mod,         only : DRIFTER_U_CURRENT_COMPONENT, DRIFTER_TEMPERATURE,       &
                                 DRIFTER_V_CURRENT_COMPONENT, FLOAT_TEMPERATURE,         &
                                 FLOAT_SALINITY
use obs_def_mod,          only : obs_def_type, set_obs_def_time, set_obs_def_key,        &
                                 set_obs_def_error_variance, set_obs_def_location,       &
                                 set_obs_def_type_of_obs
use parse_args_mod,       only : get_args_from_string 

implicit none

character(len=*), parameter :: source = 'bmkg_to_obs'

! Parameters
integer, parameter  :: PRODUCT_UNKNOWN = 0
integer, parameter  :: PRODUCT_SVP     = 1
integer, parameter  :: PRODUCT_ARVORC  = 2
integer, parameter  :: PRODUCT_ARVORI  = 3
integer, parameter  :: MAX_NUM_FIELDS  = 1000
integer, parameter  :: MAX_FIELDS_LEN  = 15000 

character(len=*), parameter :: EMPTY_ENTRY = '_EMPTY_'

character(len=512) :: string1, string2

! File variables
character(len=256) :: next_infile
integer            :: io, iunit, filenum
integer            :: num_new_obs, nfiles
logical            :: from_list = .false., first_obs = .true.

type insitu_database
    integer   :: prod   ! product type: arvorc, arvori, svp
    integer   :: nrows  ! number of data entries (including header)
    integer   :: ncols  ! number of data fields
    integer   :: nobs   ! number of valid obs in the database
    character :: del    ! delimiter style 
    
    character(len=512) :: fields(MAX_NUM_FIELDS)  ! header field names 
end type insitu_database

type(insitu_database) :: db 


! Obs sequence
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time
integer                 :: num_copies = 1,   &   ! number of copies in sequence
                           num_qc     = 1        ! number of QC entries
real(r8)                :: obs_qc

! In-situ data types
type arvor
    ! ARVOR profile
    type(time_type) :: dat 

    real(r8) :: lat
    real(r8) :: lon
    real(r8) :: dep
    real(r8) :: temp
    real(r8) :: sal 
end type arvor

type svp
    ! SVP surface drifters
    type(time_type) :: dat

    real(r8) :: lat
    real(r8) :: lon
    real(r8) :: dep
    real(r8) :: sst
    real(r8) :: uvel
    real(r8) :: vvel
end type svp

type(arvor), allocatable :: argo(:)
type(svp),   allocatable :: drift(:)

integer,  dimension(2)   :: float_types   = [FLOAT_TEMPERATURE, FLOAT_SALINITY]
integer,  dimension(3)   :: drifter_types = [DRIFTER_TEMPERATURE, DRIFTER_U_CURRENT_COMPONENT, &
                                             DRIFTER_V_CURRENT_COMPONENT] 
real(r8), dimension(2)   :: float_errs  
real(r8), dimension(3)   :: drifter_errs

!------------------------------------------------------------------------
!  Declare namelist parameters
character(len=256) :: file_in           = ''
character(len=256) :: file_list         = ''
character(len=256) :: file_out          = 'obs_seq.bmkg'
integer            :: avg_obs_per_file  = 500000
logical            :: convert_svp       = .true. 
logical            :: convert_arvorc    = .true. 
logical            :: convert_arvori    = .true.
real(r8)           :: obs_error_tmp     = 0.02_r8  ! Profiling floats are highly accurate
real(r8)           :: obs_error_sal     = 0.02_r8  ! ARVOR CTDs are precise
real(r8)           :: obs_error_vel     = 0.10_r8  ! Can vary locally
real(r8)           :: obs_error_sst     = 0.20_r8  ! Derived from small sensors, slightly higher noise
logical            :: debug             = .true.


namelist /bmkg_to_obs_nml/ file_in,          &
                           file_list,        &
                           convert_svp,      &
                           convert_arvorc,   &
                           convert_arvori,   &
                           obs_error_tmp,    &
                           obs_error_sal,    &
                           obs_error_sst,    &
                           obs_error_vel,    &
                           file_out,         &
                           avg_obs_per_file, &
                           debug

! Start Converter
call initialize_utilities()

! Read the namelist options
call find_namelist_in_file('input.nml', 'bmkg_to_obs_nml', iunit)
read(iunit, nml = bmkg_to_obs_nml, iostat = io)

if (do_nml_file()) write(nmlfileunit, nml=bmkg_to_obs_nml)
if (do_nml_term()) write(     *     , nml=bmkg_to_obs_nml)

! Set the calendar kind
call set_calendar_type(GREGORIAN)

! Check the files
if (file_in   /= '' .and. file_list /= '') then
   string1 = 'One of input HF file or the file list must be NULL'
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
float_errs    = [obs_error_tmp, obs_error_sal]
drifter_errs  = [obs_error_sst, obs_error_vel, obs_error_vel]

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
call set_copy_meta_data(obs_seq, num_copies, 'In-situ observation')
call set_qc_meta_data(obs_seq, num_qc, 'In-situ QC')

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

  ! Open it and find the product at hand
  call detect_product(next_infile, db) 

  if (.not. product_enabled(db%prod)) then
     if (debug) write(*,'(A)') 'Skipping (disabled in namelist).'
     cycle FILELOOP
  endif

  ! Start parsing
  select case (db%prod)
   case (PRODUCT_SVP)
      allocate(drift(db%nrows-1))
      call parse_svp(next_infile, db, drift)
      call add_svp(obs, prev_obs, prev_time, obs_qc, db, drift)

   case (PRODUCT_ARVORC, PRODUCT_ARVORI)
      allocate(argo(db%nrows-1))
      call parse_arvor_profile(next_infile, db, argo)
      call add_arvor_profile(obs, prev_obs, prev_time, obs_qc, db, argo)

   case default
      if (debug) write(*,'(A)') 'Unknown product; skipping.'
  endselect
  call cleanup

  filenum = filenum + 1
enddo FILELOOP

call finish_obs()
call finalize_utilities()

contains


!----------------------------------------------------------------
! Given surface drifter obs, add the metadata in the obs sequence  
subroutine add_svp(obs, prev_obs, prev_time, oqc, catalog, surf)

type(obs_type),        intent(inout) :: obs, prev_obs
type(time_type),       intent(inout) :: prev_time
real(r8),              intent(in)    :: oqc
type(insitu_database), intent(in)    :: catalog
type(svp),             intent(in)    :: surf(:)

type(time_type)    :: odat
real(r8)           :: olon, olat, odep, oval, oerr
integer            :: otype, iobs, osec, oday, k

do k = 1, 3
   otype = drifter_types(k)
   oerr  = drifter_errs(k)  
 
   do iobs = 1, catalog%nobs
      olon = surf(iobs)%lon
      olat = surf(iobs)%lat
      odep = surf(iobs)%dep
      odat = surf(iobs)%dat

      if (k == 1) oval = surf(iobs)%sst
      if (k == 2) oval = surf(iobs)%uvel
      if (k == 3) oval = surf(iobs)%vvel

      call get_time(odat, osec, oday)
      call create_3d_obs(olat, olon, odep, VERTISSURFACE, oval, otype, oerr, oday, osec, oqc, obs)
      call add_obs_to_seq(obs_seq, obs, odat, prev_obs, prev_time, first_obs)
   enddo
   ! OBS
enddo

end subroutine add_svp


!------------------------------------------------------
! Given a profile, add the metadata in the obs sequence  
subroutine add_arvor_profile(obs, prev_obs, prev_time, oqc, catalog, prof)

type(obs_type),        intent(inout) :: obs, prev_obs
type(time_type),       intent(inout) :: prev_time
real(r8),              intent(in)    :: oqc
type(insitu_database), intent(in)    :: catalog
type(arvor),           intent(in)    :: prof(:)

type(time_type)    :: odat
real(r8)           :: olon, olat, odep, oval, oerr
integer            :: otype, iobs, osec, oday, k

do k = 1, 2
   otype = float_types(k)
   oerr  = float_errs(k)

   do iobs = 1, catalog%nobs
      olon = prof(iobs)%lon
      olat = prof(iobs)%lat
      odep = prof(iobs)%dep
      odat = prof(iobs)%dat

      if (k == 1) oval = prof(iobs)%temp
      if (k == 2) oval = prof(iobs)%sal

      call get_time(odat, osec, oday)
      call create_3d_obs(olat, olon, odep, VERTISHEIGHT, oval, otype, oerr, oday, osec, oqc, obs)
      call add_obs_to_seq(obs_seq, obs, odat, prev_obs, prev_time, first_obs)
   enddo
   ! OBS
enddo

end subroutine 


!---------------------------------------------------
! Read the data profile and the associated dates
subroutine parse_svp(fname, catalog, surf)

character(len=*), parameter :: routine = 'parse_svp'

character(*),          intent(in)    :: fname
type(insitu_database), intent(inout) :: catalog
type(svp),             intent(out)   :: surf(:)

integer              :: i, k, nfields
integer              :: idat, ilat, ilon, isst, iuvel, ivvel
real(r8)             :: lat, lon, sst, uvel, vvel
character(len=10000) :: line
character(len=512)   :: entries(MAX_NUM_FIELDS)

iunit = open_file(fname, action='read', form='formatted')

! Skip the header
read(iunit, '(A)', iostat=io) line
if (io /= 0) then
   if (debug) print *, 'Got bad read code from input file, io = ', io
   call close_file(iunit)
   return
endif

! Indices of the fields of interest
idat  = find_field(catalog%ncols, catalog%fields, 'dat')
ilat  = find_field(catalog%ncols, catalog%fields, 'lat')
ilon  = find_field(catalog%ncols, catalog%fields, 'lon')
isst  = find_field(catalog%ncols, catalog%fields, 'sea_temperature')
iuvel = find_field(catalog%ncols, catalog%fields, 'u')
ivvel = find_field(catalog%ncols, catalog%fields, 'v')

! Read actual data
k = 0

OBSLOOP: do i = 1, catalog%nrows-1
   read(iunit, '(A)') line

   ! start reading each row of data points
   call split_fields(line, catalog%del, nfields, entries)

   ! just in case we get a bad line
   if (minval([idat, ilat, ilon, isst, iuvel, ivvel]) < 1 .or. &
       maxval([idat, ilat, ilon, isst, iuvel, ivvel]) > nfields) cycle OBSLOOP

   ! location and physical data
   lat  = string_to_real(entries(ilat))
   lon  = string_to_real(entries(ilon))
   sst  = string_to_real(entries(isst))
   uvel = string_to_real(entries(iuvel))
   vvel = string_to_real(entries(ivvel))

   if (lat  == MISSING_R8 .or. &
       lon  == MISSING_R8 .or. &
       sst  == MISSING_R8 .or. &
       uvel == MISSING_R8 .or. &
       vvel == MISSING_R8 ) cycle OBSLOOP

   k = k+1
   if (lon < 0.0_r8) lon = lon + 360.0_r8

   ! time
   call parse_time(entries(idat), surf(k)%dat)

   surf(k)%lon  = lon
   surf(k)%lat  = lat
   surf(k)%dep  = 0.0_r8          ! surface
   surf(k)%sst  = sst - t_kelvin  ! celsius
   surf(k)%uvel = uvel
   surf(k)%vvel = vvel

   if (debug) then
      write(string1, '(A, i3, 6(A,f10.4), A)') ' * obs #',   k , ', lat:  ', surf(k)%lat, &
                                       ', lon:'  , surf(k)%lon , ', dep:  ', surf(k)%dep,  &
                                       ', SST:'  , surf(k)%sst , ', Uvel: ', surf(k)%uvel, &
                                       ', Vvel: ', surf(k)%vvel, ', date: '

      call print_date(surf(k)%dat, str=string1)
   endif
enddo OBSLOOP

! Enter the number of acceptable observations
catalog%nobs = k

call close_file(iunit)

end subroutine parse_svp


!---------------------------------------------------
! Read the data profile and the associated dates
subroutine parse_arvor_profile(fname, catalog, prof)

character(len=*), parameter :: routine = 'parse_arvor_profile'

character(*),          intent(in)    :: fname
type(insitu_database), intent(inout) :: catalog 
type(arvor),           intent(out)   :: prof(:)

integer              :: i, k, nfields
integer              :: idat, ilat, ilon, ipres, itemp, isaln
real(r8)             :: lat, lon, pres, temp, sal
character(len=10000) :: line
character(len=512)   :: entries(MAX_NUM_FIELDS)

iunit = open_file(fname, action='read', form='formatted')

! Skip the header
read(iunit, '(A)', iostat=io) line
if (io /= 0) then
   if (debug) print *, 'Got bad read code from input file, io = ', io
   call close_file(iunit)
   return
endif

! Indices of the fields of interest
idat  = find_field(catalog%ncols, catalog%fields, 'dat')
ilat  = find_field(catalog%ncols, catalog%fields, 'lat')
ilon  = find_field(catalog%ncols, catalog%fields, 'lon')
ipres = find_field(catalog%ncols, catalog%fields, 'sea_water_pressure')
itemp = find_field(catalog%ncols, catalog%fields, 'sea_water_temperature')
isaln = find_field(catalog%ncols, catalog%fields, 'salinity') 

! Read actual data
k = 0

OBSLOOP: do i = 1, catalog%nrows-1 
   read(iunit, '(A)') line

   ! start reading each row of data points
   call split_fields(line, catalog%del, nfields, entries)

   if (minval([idat, ilat, ilon, ipres, itemp, isaln]) < 1 .or. &
       maxval([idat, ilat, ilon, ipres, itemp, isaln]) > nfields) cycle OBSLOOP
 
   ! location and physical data
   lat  = string_to_real(entries(ilat))
   lon  = string_to_real(entries(ilon))
   pres = string_to_real(entries(ipres))
   temp = string_to_real(entries(itemp))
   sal  = string_to_real(entries(isaln))
    
   if (lat  == MISSING_R8 .or. &
       lon  == MISSING_R8 .or. &
       pres == MISSING_R8 .or. &
       temp == MISSING_R8 .or. & 
       sal  == MISSING_R8 ) cycle OBSLOOP

   k = k+1  
   if (lon < 0.0_r8) lon = lon + 360.0_r8 

   ! time
   call parse_time(entries(idat), prof(k)%dat)
    
   prof(k)%lon  = lon 
   prof(k)%lat  = lat
   prof(k)%dep  = depth_from_pressure(pres, lat)
   prof(k)%temp = temp - t_kelvin
   prof(k)%sal  = sal 

   if (debug) then
      write(string1, '(A, i3, 5(A,f10.4), A)') ' * obs #', k,   ', lat:', prof(k)%lat, &
                                    ', lon:', prof(k)%lon,  ', dep:', prof(k)%dep, &
                                    ', T:'  , prof(k)%temp, ', S:'  , prof(k)%sal, &
                                    ', date: '
      
      call print_date(prof(k)%dat, str=string1)
   endif  
enddo OBSLOOP

! Enter the number of acceptable observations
catalog%nobs = k

call close_file(iunit)

end subroutine parse_arvor_profile


!------------------------------------------------
! Adapt get_args_from_string()
subroutine split_fields(line, delim, nfields, fields)
    
character(len=*), intent(in)  :: line
character,        intent(in)  :: delim
integer,          intent(out) :: nfields
character(len=*), intent(out) :: fields(:)

character(len=MAX_FIELDS_LEN) :: work

! Clean the line then parse it through 
! one of our utilities.
work = normalize_delims(line, delim)
call get_args_from_string(work, nfields, fields)

end subroutine split_fields


!----------------------------------------------------------------------
! Replace ',' and ';' with blanks to use internal parser
! We also need to treat empty fields so that we don't
! collapse with the spaces and cause any column drifts. 
function normalize_delims(line, delim) result(out_line)
 
character(len=*), intent(in)  :: line
character,        intent(in)  :: delim

character(len=MAX_FIELDS_LEN) :: out_line
integer                       :: i, j, L, k, lee
logical                       :: prev_is_delim

! Start as with a delimiter 
out_line      = ' '
prev_is_delim = .true.
    
j = 1
L = len_trim(line)

lee = len(EMPTY_ENTRY)

! Go over the line 1 character at a time
do i = 1, L
   if (line(i:i) == char(13)) cycle
   if (line(i:i) == delim) then
      ! Found a delim
      if (prev_is_delim) then
         ! insert placeholder + 1 space
         out_line(j:j+lee-1) = EMPTY_ENTRY
         j = j+lee
         out_line(j:j) = ' '
         
         j = j+1
      else
         ! normal delimiter
         out_line(j:j) = ' '
         j = j+1
      endif
      prev_is_delim = .true.
      if (j > MAX_FIELDS_LEN - 64) exit ! prevent overflow; 64 is a small cushion
   else
      out_line(j:j) = line(i:i) 

      j = j+1
      prev_is_delim = .false.
      if (j > MAX_FIELDS_LEN - 64) exit
   endif
enddo

! Trailing empty field: line ends with a delimiter (or several)
if (L > 0 .and. line(L:L) == delim) then
   out_line(j:j+lee-1) = EMPTY_ENTRY
   j = j + lee
endif

! Trim right spaces
k = j - 1
do while (k >= 1 .and. out_line(k:k) == ' ')
   k = k - 1
enddo

if (k < 1) then
   out_line = ''
else
   out_line = out_line(1:k)
endif
 
end function normalize_delims


!------------------------------------------------------------
! Get index of requested field
integer function find_field(nfields, fields, key) result(idx)

integer           :: nfields
character(len=*)  :: fields(:)
character(len=*)  :: key

integer           :: i
character(len=64) :: field_name, field_key

field_key = adjustl(trim(key))
call to_upper(field_key)

idx = -1

do i = 1, nfields
   field_name = adjustl(trim(fields(i)))
   call to_upper(field_name)
   if (field_name == field_key) then
      idx = i
      return
   endif
enddo

end function find_field


!-----------------------------------------
! Using the date string, set the DART time
subroutine parse_time(datestr, obs_time)

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

end subroutine parse_time


!-------------------------------------------------------
! Lightweight in-situ platform detector before parsing.
subroutine detect_product(fname, catalog) 

character(len=*), parameter :: routine = 'detect_product'

character(*),          intent(in)  :: fname
type(insitu_database), intent(out) :: catalog

integer              :: ip, nfields
character(len=10000) :: line
character(len=512)   :: entries(MAX_NUM_FIELDS)

! Set it to an unknown float/drifter
catalog%prod = PRODUCT_UNKNOWN

! Number of data entries in file
call find_textfile_dims(fname, catalog%nrows)

iunit = open_file(fname, action='read', form='formatted')

! Read header
read(iunit, '(A)', iostat=io) line
if(io /= 0) then
   if (debug) print *, 'Got bad read code from input file, io = ', io
   call close_file(iunit)
   return
endif
catalog%del = detect_delim(line)  ! is it a ',' or a ';'

! Split header fields
call split_fields(line, catalog%del, catalog%ncols, catalog%fields)

! Get 'process' position
ip = find_field(catalog%ncols, catalog%fields, 'process')
if (ip < 1) then
   call close_file(iunit)
   catalog%prod = PRODUCT_UNKNOWN
   return
endif

! Obtain the product type
read(iunit, '(A)', iostat=io) line
if(io /= 0) then
   if (debug) print *, 'Got bad read code from input file, io = ', io
   call close_file(iunit)
   return
endif

call split_fields(line, catalog%del, nfields, entries)
catalog%prod = classify_insitu_data(entries(ip))

call close_file(iunit)

end subroutine detect_product


!-----------------------------------------------------
! Obtain the data platform from the process name
integer function classify_insitu_data(field) result(platform)

character(*) :: field

field = adjustl(trim(field))
call to_upper(field)

if (field == 'ARVORC') then 
   platform = PRODUCT_ARVORC 
elseif (field == 'ARVORI') then 
   platform = PRODUCT_ARVORI 
elseif (field == 'DRIFTER' .or. field == 'SVP') then 
   platform = PRODUCT_SVP
else
   platform = PRODUCT_UNKNOWN
endif

end function classify_insitu_data


!--------------------------------------
! Check whether a product is enabled
logical function product_enabled(prod)

integer, intent(in) :: prod

select case (prod)
  case (PRODUCT_SVP)
     product_enabled = convert_svp
  case (PRODUCT_ARVORC)
     product_enabled = convert_arvorc
  case (PRODUCT_ARVORI)
     product_enabled = convert_arvori
  case default
     product_enabled = .false.
endselect

end function product_enabled


!-------------------------------------------------
! Because it's a CSV file, we should expect a ','
! but the sample data file came with a ';' so
! let's take care of both cases
function detect_delim(line) result(delim)

character(len=*), intent(in) :: line

character :: delim
integer   :: i, nsemi, ncomma, L

L = len_trim(line)

nsemi  = 0
ncomma = 0

do i = 1, L
   if (line(i:i) == ';') nsemi  = nsemi  + 1
   if (line(i:i) == ',') ncomma = ncomma + 1
enddo
  
! Is there a better check? 
if (nsemi >= ncomma) then
   delim = ';'
else
  delim = ','
endif

end function detect_delim


!------------------------------------------------------------
! All collected obs (if any) are written in the seq file. 
subroutine finish_obs()

integer :: obs_num

obs_num = get_num_obs(obs_seq)

if (obs_num > 0) then
   print *
   if (debug) write(*, '(A, i0, A)') '>>>> Ready to write ', obs_num, ' observations:'

   call write_obs_seq(obs_seq, file_out)
   call destroy_obs(obs)
   call destroy_obs_sequence(obs_seq)
else
   string1 = 'No obs were converted.'
   call error_handler(E_ERR, source, string1)
endif
call error_handler(E_MSG, source, 'Finished successfully.')

end subroutine finish_obs


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
g    = 9.780318_r8 * (1._r8 + 5.2788e-3_r8*sin2 + 2.36e-5_r8*sin4) &
       + 1.092e-6_r8 * p_dbar

! Polynomial (UNESCO 1983): depth (m) from pressure (dbar)
! depth = (((-1.82e-15*p + 2.279e-10)*p - 2.2512e-5)*p + 9.72659)*p / g
depth_m = (((-1.82e-15_r8*p_dbar + 2.279e-10_r8)*p_dbar - 2.2512e-5_r8)*p_dbar + 9.72659_r8) &
          * p_dbar / g
end function depth_from_pressure


!------------------------------------------------------------
! Clear up memory 
subroutine cleanup()

if (allocated(argo) ) deallocate(argo)
if (allocated(drift)) deallocate(drift)

end subroutine cleanup

end program bmkg_to_obs

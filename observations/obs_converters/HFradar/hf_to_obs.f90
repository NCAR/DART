! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


! Obs converter for High-Frequency (HF) radar data. 
! The code supports surface ocean current vectors in both eastward and 
! northward directions. It also supports radial currents (line-of-sight).

! For the gridded TOTAL (vector) surface currents: 
! - Variables: u and v are directly observed 
! - Source   : HF Radar Derived Surface Currents obtained from CODAR combine 
!              method (i.e., totals created by fusing more than 2 radials).
! - Origin   : DALI and BATU radial sites
! - Grid     : Regular (e.g., 1 km with separate 1D lat & lon)
! - Extras   : Obs errors and error covariances too

! For the single site radial currents: 
! - Variables: velo - radial sea water velocity away from the instrument)
!              head - direction away from instrument
!              bear - phi; bearing from instrument (to each cell) 
! - Source   : Near real time surface ocean radial velocity at "BATU"
! - Grid     : 2D curvilinear lat(x,y), lon(x,y) distances from instrument
! - FO       : Radial Velocity = U cos(phi) + V sin(phi) 

program hf_to_obs

use types_mod,            only : r8, i8, i4, digits12, deg2rad, rad2deg
use time_manager_mod,     only : time_type, set_calendar_type, GREGORIAN, set_time,      &
                                 get_time, print_time, set_date, get_date, print_date,   &
                                 write_time, operator(+), operator(-)
use utilities_mod,        only : initialize_utilities, find_namelist_in_file, E_ERR,     &
                                 check_namelist_read, nmlfileunit, error_handler, E_MSG, &
                                 finalize_utilities, do_nml_file, get_next_filename,     &
                                 do_nml_term, find_textfile_dims, file_exist
use location_mod,         only : VERTISSURFACE, set_location
use obs_sequence_mod,     only : obs_type, obs_sequence_type, init_obs, get_num_obs,     &
                                 static_init_obs_sequence, init_obs_sequence,            &
                                 set_copy_meta_data, set_qc_meta_data, write_obs_seq,    & 
                                 destroy_obs_sequence, insert_obs_in_seq, set_qc,        &
                                 set_obs_values, set_obs_def, destroy_obs
use obs_utilities_mod,    only : add_obs_to_seq
use obs_kind_mod,         only : HFRADAR_U_CURRENT_COMPONENT, HFRADAR_RADIAL_VELOCITY,   &
                                 HFRADAR_V_CURRENT_COMPONENT
use obs_def_mod,          only : obs_def_type, set_obs_def_time, set_obs_def_key,        &
                                 set_obs_def_error_variance, set_obs_def_location,       &
                                 set_obs_def_type_of_obs
use netcdf_utilities_mod, only : nc_check, nc_open_file_readonly, nc_close_file,         &
                                 nc_get_variable, nc_get_attribute_from_variable,        &
                                 nc_get_dimension_size, nc_get_variable_size,            &
                                 nc_variable_exists, nc_get_global_attribute,            &
                                 nc_global_attribute_exists
use obs_def_ocean_mod,    only : set_hf_radial_vel
use netcdf

implicit none

character(len=*), parameter :: source = 'hf_to_obs'

integer, parameter  :: QC_FLAG_MIN = 0       ! Radial velocity: acceptable (no flags raised)
integer, parameter  :: QC_FLAG_MAX = 1024    ! Radial velocity 'beyond range' flag

real(r8), parameter :: OBS_ERROR_SD_MIN = 1.0e-3_r8  ! Lower bound for obs_error_sd
real(r8), parameter :: OBS_ERROR_SD_MAX = 4.0e-1_r8  ! Upper bound for obs_error_sd
real(r8), parameter :: CM2M             = 1.0e-2_r8  ! Conversion from cm to meters 

integer  :: num_copies = 1,   &   ! number of copies in sequence
            num_qc     = 1        ! number of QC entries

character(len=30), dimension(3) :: obs_names = ['Eastward Sea Water Velocity ' , &
                                                'Northward Sea Water Velocity', & 
                                                'Radial Sea Water Velocity   '   ]

character(len=256)      :: next_infile, instrument
character(len=512)      :: string1, string2, string3
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, prev_time

integer                 :: ncid, io, iunit, filenum
integer                 :: num_new_obs, nfiles, dummy
integer                 :: hf_kind, obs_num
integer                 :: ilon, ilat, nlat, nlon
integer                 :: ix, iy, nx, ny

logical                 :: from_list = .false., first_obs = .true.
real(r8)                :: obs_qc, rmin_km, rmax_km 
real(r8)                :: missing_u, missing_v

! Debugging counters
integer :: totals_kept     = 0
integer :: radials_kept    = 0
integer :: drop_flag       = 0       ! radials-only
integer :: drop_range      = 0       ! radials-only
integer :: drop_nan_radial = 0       ! radials-only
integer :: drop_nan_u      = 0       ! totals U
integer :: drop_nan_v      = 0       ! totals V

! HF radar sites and instrument info
integer, parameter :: MAX_SITES = 20
character(len=64)  :: site_names(MAX_SITES) = ''
integer            :: num_sites = 0

! Arrays for TOTALS; vectored U & V
real(r8), allocatable :: lon(:), lat(:)         ! Location 
real(r8), allocatable :: u(:, :), v(:, :)       ! Totals for U and V
real(r8), allocatable :: sd_u(:, :), sd_v(:, :) ! Obs error SD for U and V

! Stats for TOTALS
real(r8) :: sd_u_min = huge(1.0_r8), sd_u_max = -huge(1.0_r8), sd_u_sum = 0.0_r8
real(r8) :: sd_v_min = huge(1.0_r8), sd_v_max = -huge(1.0_r8), sd_v_sum = 0.0_r8
real(r8) :: u_min    = huge(1.0_r8), u_max    = -huge(1.0_r8), u_sum    = 0.0_r8
real(r8) :: v_min    = huge(1.0_r8), v_max    = -huge(1.0_r8), v_sum    = 0.0_r8
integer  :: num_sd_u = 0           , num_sd_v = 0
integer  :: n_uval   = 0           , n_vval   = 0

! Arrays for RADIAL velocities
real(r8), allocatable :: rlon(:, :), rlat(:, :) ! Cell geolocation
real(r8), allocatable :: rvel(:, :)             ! Measured radial velocity (away from radar) 
real(r8), allocatable :: rflg(:, :)             ! Vector (quality) flag: 0 <= good < 1024 
real(r8), allocatable :: phi(:, :)              ! Bearing (azimuth) of beam; clockwise from true north
real(r8), allocatable :: rnge(:,:)              ! Range from the instrument to that grid cell (along beam line)
real(r8), allocatable :: rstd(:,:)              ! Obs error standard deviation

! Stats for RADIALS
real(r8) :: sd_r_min = huge(1.0_r8), sd_r_max = -huge(1.0_r8), sd_r_sum = 0.0_r8
real(r8) :: r_min    = huge(1.0_r8), r_max    = -huge(1.0_r8), r_sum    = 0.0_r8
integer  :: num_sd_r = 0           , n_rval   = 0

!------------------------------------------------------------------------
!  Declare namelist parameters
character(len=256) :: file_in           = ''
character(len=256) :: file_list         = ''
character(len=256) :: file_out          = 'obs_seq.hf'
integer            :: avg_obs_per_file  = 500000
real(r8)           :: radial_rmin_km    = -1.0_r8   ! < 0 derive from data 
real(r8)           :: radial_rmax_km    = -1.0_r8   ! < 0 derive from data
logical            :: debug             = .true.


namelist /hf_to_obs_nml/ file_in,          &
                         file_list,        &   
                         file_out,         &
                         radial_rmin_km,   &
                         radial_rmax_km,   &
                         avg_obs_per_file, &
                         debug   

! Start Converter
call initialize_utilities()

! Read the namelist options
call find_namelist_in_file('input.nml', 'hf_to_obs_nml', iunit)
read(iunit, nml = hf_to_obs_nml, iostat = io) 

if (do_nml_file()) write(nmlfileunit, nml=hf_to_obs_nml)
if (do_nml_term()) write(     *     , nml=hf_to_obs_nml)

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
   call find_textfile_dims(file_list, nfiles, dummy)
   num_new_obs = avg_obs_per_file * nfiles 
endif

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
call set_copy_meta_data(obs_seq, num_copies, 'HF radar observation')
call set_qc_meta_data(obs_seq, num_qc, 'HF radar QC')

! Loop over the obs files
filenum = 1
obs_num = 1
obs_qc  = 0.0_r8

FILELOOP: do 
  ! Get the HF data file
  if (from_list) then 
     next_infile = get_next_filename(file_list, filenum) 
  else
     next_infile = file_in
     if (filenum > 1) next_infile = ''
  endif
  if (next_infile == '') exit FILELOOP  
  print * 
  if (debug) write(*, '(A, i0, X, A)') 'Input file: #', filenum, next_infile

  ! Open it and configure it
  ! hf_kind = 0: Totals, 1: Radials 
  call configure_HF_file(next_infile, ncid, hf_kind, obs_time)  

  if (hf_kind == 0) then
     ! More straight-forward case:
     ! We have access to the total U & V. 
     ! Below, we read the obs locations, currents values
     ! and associated errors 
     call read_hf_totals(ncid)
     do ilon = 1, nlon
        do ilat = 1, nlat
           ! Adding a U obs
           if (good_total_obs(u(ilon,ilat), 'U')) then
              call fill_obs(obs, prev_obs, prev_time, obs_num, lon(ilon), &
                   lat(ilat), HFRADAR_U_CURRENT_COMPONENT, obs_time,      & 
                   sd_u(ilon, ilat), u(ilon, ilat), obs_qc)
              totals_kept = totals_kept + 1
           endif
           ! Adding a V obs
           if (good_total_obs(v(ilon,ilat), 'V')) then
              call fill_obs(obs, prev_obs, prev_time, obs_num, lon(ilon), &
                   lat(ilat), HFRADAR_V_CURRENT_COMPONENT, obs_time,      &
                   sd_v(ilon, ilat), v(ilon, ilat), obs_qc)
              totals_kept = totals_kept + 1
           endif
           call bump_val_total(u(ilon,ilat), v(ilon,ilat))
           call bump_sd_total(sd_u(ilon,ilat), sd_v(ilon,ilat))  
        enddo
     enddo     
  else
     ! Need to read radial velocities and extra metadata (e.g., angle)
     ! for forward operator calculations
     call read_hf_radials(ncid)
     do ix = 1, nx
        do iy = 1, ny
           ! Adding radial velocity obs
           if (good_radial_obs(rvel(iy, ix), rflg(iy, ix), rnge(iy, ix))) then 
              call fill_obs(obs, prev_obs, prev_time, obs_num, rlon(iy, ix), & 
                   rlat(iy, ix), HFRADAR_RADIAL_VELOCITY, obs_time,          &
                   rstd(iy, ix), rvel(iy, ix), obs_qc, phi(iy, ix))
              call bump_radials(rvel(iy, ix), rstd(iy, ix))
              radials_kept = radials_kept + 1
           endif 
        enddo
     enddo
  endif 
  call nc_close_file(ncid, source)
  call cleanup()

  filenum = filenum + 1
enddo FILELOOP 

! Done with main loop; write [1] summary and [2] the obs_seq out
if (debug) call report_summary()

call finish_obs()
call finalize_utilities()


contains


!------------------------------------------------------------
! Fill in an obs: loc in 3D in this case (at the surface).
! Once filled, add it to the sequence.
subroutine fill_obs(obs, prev_obs, prev_time, onum, olon, olat, otype, otime, oerr, oval, oqc, angle)

type(obs_type),     intent(inout) :: obs, prev_obs
type(time_type),    intent(inout) :: prev_time
integer,            intent(inout) :: onum
real(r8),           intent(in)    :: olon, olat, oerr, oval, oqc
integer,            intent(in)    :: otype
type(time_type),    intent(in)    :: otime
real(r8), optional, intent(in)    :: angle

type(obs_def_type) :: obs_def
real(r8)           :: ovalue(1), oqcval(1)
integer            :: instrument_ID
integer            :: key 

ovalue(1) = oval
oqcval(1) = oqc

call set_obs_def_location(obs_def, set_location(olon, olat, 0.0_r8, VERTISSURFACE))
call set_obs_def_type_of_obs(obs_def, otype)
call set_obs_def_time(obs_def, otime)
call set_obs_def_error_variance(obs_def, oerr * oerr)

if (present(angle)) then 
   instrument_ID = string_hash_id(instrument)
   call set_hf_radial_vel(key, instrument_ID, angle)
else
   key = onum
endif

call set_obs_def_key(obs_def, key)
call set_obs_def(obs, obs_def)
call set_obs_values(obs, ovalue)
call set_qc(obs, oqcval)

onum = onum + 1  

! Add the obs to the sequence
! and handle time properly
call add_obs_to_seq(obs_seq, obs, otime, prev_obs, prev_time, first_obs)

end subroutine fill_obs


!------------------------------------------------------------
! All collected obs (if any) are written in the seq file. 
subroutine finish_obs()

if (obs_num > 0) then
   print *
   if (debug) write(*, '(A, i0, A)') '>>>> Ready to write ', get_num_obs(obs_seq), ' observations:'
   if (get_num_obs(obs_seq) > 0)  call write_obs_seq(obs_seq, file_out)

   call destroy_obs(obs)
   if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
else
   string1 = 'No obs were converted.'
   call error_handler(E_MSG, source, string1)
endif
call error_handler(E_MSG, source, 'Finished successfully.')

end subroutine finish_obs


!------------------------------------------------------------
! Read the obs location plus vectors: U and V
! in addition to the obs error std
subroutine read_hf_totals(ncid)

character(len=*), parameter :: routine = 'read_hf_totals'

integer, intent(in) :: ncid
real(r8)            :: uval, vval
character(len=64)   :: vel_unit

if (debug) call print_date(obs_time, str= '   * Reading in TOTAL U and V data:')

nlat = nc_get_dimension_size(ncid, 'lat', routine)
nlon = nc_get_dimension_size(ncid, 'lon', routine)

if (debug) then
   if (nc_get_dimension_size(ncid,'time', routine) /= 1) then
      string1 = '"time" dimension > 1 not supported.'
      call error_handler(E_ERR, routine, string1, source)
   endif
endif

allocate(lat(nlat), lon(nlon))
allocate(u(nlon, nlat), v(nlon, nlat))
allocate(sd_u(nlon, nlat), sd_v(nlon, nlat))

call nc_get_variable(ncid, 'lat', lat, routine)
call nc_get_variable(ncid, 'lon', lon, routine)
where(lon < 0.0_r8) lon = lon + 360.0_r8

call nc_get_attribute_from_variable(ncid, 'u', '_FillValue', missing_u, routine)
call nc_get_attribute_from_variable(ncid, 'v', '_FillValue', missing_v, routine)

! Read U, V and provided obs errors
call nc_get_variable(ncid, 'u', u, routine)
call nc_get_variable(ncid, 'v', v, routine)

! Make sure the errors are in the netcdf file
sd_u = OBS_ERROR_SD_MIN
if (nc_variable_exists(ncid,'stdu')) then
   call nc_get_variable(ncid, 'stdu', sd_u,routine)
endif

sd_v = OBS_ERROR_SD_MIN
if (nc_variable_exists(ncid,'stdv')) then
   call nc_get_variable(ncid,'stdv', sd_v, routine)
endif

if (debug) write(*, '(3X, A)') '* Read in U, V, stdU, and stdV from the input obs file.'

! Make sure we get u and v in m/s
call nc_get_attribute_from_variable(ncid, 'u', 'units', vel_unit)
call velocity_units(adjustl(vel_unit), u, obs_names(1))
call nc_get_attribute_from_variable(ncid, 'v', 'units', vel_unit)
call velocity_units(adjustl(vel_unit), v, obs_names(2))

! Adjust obs errors
do ilon = 1, nlon
   do ilat = 1, nlat
      if (sd_u(ilon,ilat) /= sd_u(ilon,ilat)) sd_u(ilon,ilat) = OBS_ERROR_SD_MIN
      if (sd_v(ilon,ilat) /= sd_v(ilon,ilat)) sd_v(ilon,ilat) = OBS_ERROR_SD_MIN
      
      uval = abs(u(ilon, ilat))
      vval = abs(v(ilon, ilat)) 
      if (uval /= missing_u .and. uval == uval) &
         sd_u(ilon, ilat) = obs_err_sd_check(sd_u(ilon, ilat))
      if (vval /= missing_v .and. vval == vval) &
         sd_v(ilon, ilat) = obs_err_sd_check(sd_v(ilon, ilat))
   enddo
enddo 

if (debug) write(*,'(3X,A)') '* Clamped obs errors to configured min/max.'

end subroutine read_hf_totals


!---------------------------------------------------------------
! Real the lin-of-sight velocity component and assign obs errors. 
subroutine read_hf_radials(ncid)

character(len=*), parameter :: routine = 'read_hf_radials'

integer, intent(in) :: ncid
real(r8)            :: sigma, r
real(r8)            :: sigma0, grow
character(len=64)   :: vel_unit, angle_unit
logical             :: attin

! Coastal CODAR radials 0.10 - 0.25
! Source: IOOS, SEACOOS guidelines
sigma0 = 0.15_r8  ! m/s, near-site error
grow   = 1.0_r8   ! scaling factor

if (debug) call print_date(obs_time, str= '   * Reading in Radial Velocity:')

nx = nc_get_dimension_size(ncid, 'x', routine)
ny = nc_get_dimension_size(ncid, 'y', routine)

allocate(rlat(ny, nx), rlon(ny, nx))
allocate(rvel(ny, nx), rflg(ny, nx))
allocate(rnge(ny, nx), phi(ny, nx))
allocate(rstd(ny, nx))

call nc_get_variable(ncid, 'lat', rlat, routine)
call nc_get_variable(ncid, 'lon', rlon, routine)
where(rlon < 0.0_r8) rlon = rlon + 360.0_r8

! Read quality flags, bearing angle, and radial velocity
call nc_get_variable(ncid, 'vflg', rflg, routine)
call nc_get_variable(ncid, 'bear',  phi, routine)
call nc_get_variable(ncid, 'velo', rvel, routine)
call nc_get_variable(ncid, 'rnge', rnge, routine)

if (debug) then
      if ( .not. ( all(shape(rlat) == shape(rlon)) .and. &
                   all(shape(rvel) == shape(rlon)) .and. &
                   all(shape(phi ) == shape(rlon)) .and. &
                   all(shape(rnge) == shape(rlon)) .and. &
                   all(shape(rflg) == shape(rlon)) ) ) then 
      string1 = 'Inconsistent 2-D dimensions'
      write(string2,'(A,2I6,A,2I6)') 'Mismatch: rlat ',shape(rlat), ' vs rlon ',shape(rlon)
      call error_handler(E_ERR, routine, string1, source, text2=string2)
   endif
endif

! Check the units of the angle from the instrument
call nc_get_attribute_from_variable(ncid, 'bear', 'units', angle_unit)
if (angle_unit(1:7) == 'degrees') then 
   ! Do nothing because the forward operator will do the conversion
   if (debug) write(*, '(5X, A)') 'Bear angle is in degrees. The FO will convert it to radians.'
elseif  (angle_unit(1:7) == 'radians') then
   phi = phi * rad2deg
   if (debug) write(*, '(5X, A)') 'Incoming bearing is in radians.'
else
   write(*, string1) 'Unknown angle units: '//trim(adjustl(angle_unit))
   call error_handler(E_ERR, routine, string1, source)
endif

! Make sure we get line-of-sight vel in m/s
call nc_get_attribute_from_variable(ncid, 'velo', 'units', vel_unit)
call velocity_units(adjustl(vel_unit), rvel, obs_names(3))

! Derive safe range limits from the data 
call range_limits()

! Compute range-dependent obs error sd
rstd = sigma0 

do ix = 1, nx
   do iy = 1, ny
      sigma = sigma0
      r     = rnge(iy, ix) ! distance in km 
      if (r == r .and. rmax_km > 0.0_r8) then 
         sigma = sigma0 * (1.0_r8 + grow * (r / rmax_km))
      endif          
      rstd(iy, ix) = obs_err_sd_check(sigma)
   enddo
enddo

if (debug) write(*, '(3X, A)') '* Computed obs errors using a range-dependent strategy.'

! Assign an instrument
attin = nc_global_attribute_exists(ncid, 'Site') 

if (attin) then 
   call nc_get_global_attribute(ncid, 'Site', instrument)
   call clarify_instrument(instrument)
else
   instrument = 'UNKNOWN'
endif
call register_instrument(instrument)

end subroutine read_hf_radials


!------------------------------------------------------------
! Figure out the type of HF data that we have 
! and parse the time of the obs 
subroutine configure_HF_file(hf_file, id, hf_kind, obs_time)

character(len=*), parameter  :: routine = 'configure_HF_file'

character(len=*), intent(in) :: hf_file
integer, intent(out)         :: id 
integer, intent(out)         :: hf_kind
type(time_type), intent(out) :: obs_time

integer           :: year, month, day, hour, minute, second, ios
real(digits12)    :: time_since_init
character(len=64) :: datestr
type(time_type)   :: base_date

integer(i8) :: big_integer
integer     :: some_seconds, some_days

id = nc_open_file_readonly(hf_file)
hf_kind = 0 ! HF U, V Totals

! Check if 'Bearing From Instrument' is available
if (nc_variable_exists(id, 'bear')) hf_kind = 1   ! HF Radials case 

! Manage time
call nc_get_variable(id, 'time', time_since_init)
call nc_get_attribute_from_variable(id, 'time', 'units', datestr)

! Parse time units 
if(datestr(1:13) == 'seconds since') then 
  ! double time(time) ;
  !        time:units = "seconds since 1970-01-01 00:00:00" ;

  read(datestr, '(14x, i4, 5(1x,i2))', iostat=ios) year, month, day, hour, minute, second
  if (ios /= 0) then 
     write(string1, *) 'Unable to read time variable units. Error status was ', ios
     call error_handler(E_ERR, routine, string1, source)
  endif
  big_integer = int(time_since_init, i8)

elseif (datestr(1:11) == 'hours since') then 
  ! int time(time) ;
  !     time:units = "hours since 2025-10-06 03:10:00" ;

  read(datestr, '(12x, i4, 5(1x,i2))', iostat=ios) year, month, day, hour, minute, second
  if (ios /= 0) then 
     write(string1, *) 'Unable to read time variable units. Error status was ', ios
     call error_handler(E_ERR, routine, string1, source)
  endif
  big_integer = int(time_since_init * 3600.0_digits12, i8)  ! hours -> seconds

else
  ! Unknown time format

  write(string1, *) 'expecting time attribute units of "seconds since ..." -OR-'
  write(string2, *) '                              "hours since ..."'
  write(string3, *) 'got "'//trim(datestr)//'"'
  call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
endif

base_date    = set_date(year, month, day, hour, minute, second)
some_days    = big_integer / (24*60*60)
big_integer  = big_integer - (some_days * (24*60*60))
some_seconds = int(big_integer, i4)

obs_time = base_date + set_time(some_seconds, some_days)

end subroutine configure_HF_file


!------------------------------------------------------------
! Unit conversion 
subroutine velocity_units(units, var, oname)

character(len=*), parameter  :: routine = 'velocity_units'

character(len=*), intent(in) :: units, oname
real(r8), intent(inout)      :: var(:, :)

if (index(units,'cm') > 0 .and. index(units,'s') > 0) then
  if (debug) write(*, '(5X, A)') 'Obs "'//trim(oname)//'" unit is in '//trim(units)//'; converting to m/s.'
  var = var * CM2M

elseif (index(units,'m')>0 .and. index(units,'s')>0) then
  ! ok 
  if (debug) write(*, '(5X, A)') 'Obs "'//trim(oname)//'" unit is in '//trim(units)//'; no conversion needed.'

else
  string1 = 'Unknown units for observed currents. Exiting.'
  call error_handler(E_ERR, routine, string1, source)
endif

end subroutine velocity_units


!------------------------------------------------------------
! Use data to obtain safeguard conditions on the range of the instrument
subroutine range_limits()

real(r8) :: rmin_data, rmax_data

! Data-driven finite extents
rmin_data = minval(rnge, mask=(rnge==rnge))
rmax_data = maxval(rnge, mask=(rnge==rnge))

! Fallback if all NaN
if (rmin_data /= rmin_data .or. rmax_data /= rmax_data) then
   rmin_km = 0.0_r8
   rmax_km = huge(1.0_r8)
else
   ! Defaults: skip the very-near bins a bit; allow up to max finite range
   ! Use the user-provided rmin, rmax if provided 
   ! If not, use those from the data and skip the first 'noisy' ~1km
   rmin_km = merge(radial_rmin_km, max(0.0_r8, rmin_data + 1.0_r8), radial_rmin_km >= 0.0_r8)
   rmax_km = merge(radial_rmax_km, rmax_data,                       radial_rmax_km >= 0.0_r8)
   if (rmax_km <= rmin_km) rmax_km = rmin_km + 1.0_r8
endif

if (debug) write(*,'(5x,"Using radial range window: [",f6.2," km, ",f6.2," km]")') rmin_km, rmax_km

end subroutine range_limits


!------------------------------------------------------------
! Figure out whether we can use this total obs or not
logical function good_total_obs(val, which)
   
real(r8),           intent(in) :: val
character(len=*),   intent(in) :: which    ! 'U' or 'V'

good_total_obs = .true.

if (which(1:1) == 'U') then
   if (val /= val .or. val == missing_u) then
      drop_nan_u = drop_nan_u + 1
      good_total_obs = .false.
   endif
else
   if (val /= val .or. val == missing_v) then  
      drop_nan_v = drop_nan_v + 1
      good_total_obs = .false.
   endif
endif

end function good_total_obs


!------------------------------------------------------------
! Figure out whether we can use this radial obs or not
! Also updates run-wide counters in host scope.
logical function good_radial_obs(vel, flg, rng) 

real(r8), intent(in) :: vel, flg, rng
logical              :: ok_val, ok_flag, ok_rng

ok_val = (vel == vel) .and. (flg == flg) .and. (rng == rng)
if (.not. ok_val) then
   drop_nan_radial = drop_nan_radial + 1
   good_radial_obs = .false.
   return
endif

ok_flag = (flg >= QC_FLAG_MIN) .and. (flg < QC_FLAG_MAX)
if (.not. ok_flag) then
   drop_flag = drop_flag + 1
   good_radial_obs = .false.
   return
endif

ok_rng = (rng >= rmin_km) .and. (rng <= rmax_km)
if (.not. ok_rng) then
   drop_range = drop_range + 1
   good_radial_obs = .false.
   return
endif

good_radial_obs = .true.

end function good_radial_obs


!------------------------------------------------------------
! Assign upper and lower bounds on the obs_error_sd
function obs_err_sd_check(sd_in) result(sd_out)

real(r8) :: sd_in, sd_out

sd_out = sd_in

sd_out = max(sd_out, OBS_ERROR_SD_MIN) ! Lower bound
sd_out = min(sd_out, OBS_ERROR_SD_MAX) ! Upper bound

end function obs_err_sd_check


!------------------------------------------------------------
! Keep a record of the TOTALS valuess for summary purposes
subroutine bump_val_total(u_val, v_val)

real(r8), intent(in) :: u_val, v_val

if (u_val == u_val .and. u_val /= missing_u) then
   u_min  = min(u_min, u_val)
   u_max  = max(u_max, u_val)
   u_sum  = u_sum + u_val
   n_uval = n_uval + 1
endif

if (v_val == v_val .and. v_val /= missing_v) then
   v_min  = min(v_min, v_val)
   v_max  = max(v_max, v_val)
   v_sum  = v_sum + v_val
   n_vval = n_vval + 1
endif

end subroutine bump_val_total


!------------------------------------------------------------
! Keep a record of the TOTALS errors for summary purposes
subroutine bump_sd_total(sd_u_val, sd_v_val)

real(r8), intent(in) :: sd_u_val, sd_v_val

if (sd_u_val == sd_u_val) then
   sd_u_min = min(sd_u_min, sd_u_val)
   sd_u_max = max(sd_u_max, sd_u_val)
   sd_u_sum = sd_u_sum + sd_u_val
   num_sd_u = num_sd_u + 1
endif

if (sd_v_val == sd_v_val) then
   sd_v_min = min(sd_v_min, sd_v_val)
   sd_v_max = max(sd_v_max, sd_v_val)
   sd_v_sum = sd_v_sum + sd_v_val
   num_sd_v = num_sd_v + 1
endif

end subroutine bump_sd_total


!------------------------------------------------------------
! Keep a record of the TOTALS valuess for summary purposes
subroutine bump_radials(rvel, sd_val)

real(r8), intent(in) :: rvel, sd_val

r_min  = min(r_min, rvel)
r_max  = max(r_max, rvel)
r_sum  = r_sum + rvel
n_rval = n_rval + 1

sd_r_min = min(sd_r_min, sd_val)
sd_r_max = max(sd_r_max, sd_val)
sd_r_sum = sd_r_sum + sd_val
num_sd_r = num_sd_r + 1

end subroutine bump_radials


!------------------------------------------------------------
! Hash: Get a numeric key for the instrument 
integer function string_hash_id(str)

character(len=*), intent(in) :: str
integer :: i, h

h = 0
do i = 1, len_trim(str)
   h = iachar(str(i:i)) + 31*h
enddo
string_hash_id = abs(mod(h, 10000)) + 1

end function string_hash_id


!------------------------------------------------------------
! Remove any weird quotes or slashes from the instrument name
subroutine clarify_instrument(str)

character(len=*), intent(inout) :: str
integer :: i

! Remove leading/trailing spaces and quotes
str = adjustl(str)

! Remove any embedded double quotes or backslashes
do i = 1, len_trim(str)
   if (str(i:i) == '"'  .or. &
       str(i:i) == '''' .or. &
       str(i:i) == '\') str(i:i) = ' '
enddo

str = trim(adjustl(str))

end subroutine clarify_instrument


!------------------------------------------------------------
! Keep track of all instruments
subroutine register_instrument(name)

character(len=*), intent(in) :: name
integer :: i

! Check if this site already exists
do i = 1, num_sites
   if (trim(site_names(i)) == trim(name)) return
end do

! If not found, add it
if (num_sites < MAX_SITES) then
   num_sites = num_sites + 1
   site_names(num_sites) = trim(name)
else
   write(*,'(A)') 'Warning: exceeded MAX_SITES, some instruments ignored.'
endif

end subroutine register_instrument


!------------------------------------------------------------
! Summary info and statistics
subroutine report_summary()

character(len=512) :: site_line
integer            :: i
real(r8)           :: pct_totals, pct_radials

pct_totals = 0.0_r8
if ((totals_kept + drop_nan_u + drop_nan_v) > 0) &
   pct_totals = 100.0_r8 * totals_kept / &
   real(totals_kept + drop_nan_u + drop_nan_v, r8)

pct_radials = 0.0_r8
if ((radials_kept + drop_flag + drop_range + drop_nan_radial) > 0) &
   pct_radials = 100.0_r8 * radials_kept / &
   real(radials_kept + drop_flag + drop_range + drop_nan_radial, r8)

print *

write(*,'(5X,A)') '--- HF converter summary (all files) ---'
write(*,'(9X,A,I0," (",F5.2,"%)")') 'Totals kept     : ', totals_kept, pct_totals
write(*,'(9X,A,I0," (",F5.2,"%)")') 'Radials kept    : ', radials_kept, pct_radials
write(*,'(9X,A,I0)') 'Dropped by flag : ', drop_flag
write(*,'(9X,A,I0)') 'Dropped by range: ', drop_range
write(*,'(9X,A,I0)') 'NaN drops (U)   : ', drop_nan_u
write(*,'(9X,A,I0)') 'NaN drops (V)   : ', drop_nan_v
write(*,'(9X,A,I0)') 'NaN drops (rad) : ', drop_nan_radial

write(*, '(9X, A)') '+++++++++++++++++'
if (n_uval > 0) then
   write(*,'(9X,A,3ES12.4)') '   U   value min/max/mean (m/s)   : ', u_min, u_max, u_sum/n_uval
endif
if (n_vval > 0) then
   write(*,'(9X,A,3ES12.4)') '   V   value min/max/mean (m/s)   : ', v_min, v_max, v_sum/n_vval
endif
if (n_rval > 0) then
   write(*,'(9X,A,3ES12.4)') 'Radial value min/max/mean (m/s)   : ', r_min, r_max, r_sum/n_rval
endif

write(*, '(9X, A)') '+++++++++++++++++'
if (num_sd_u > 0) then
   write(*,'(9X,A,3ES12.4)') '   U   error SD min/max/mean (m/s): ', sd_u_min, sd_u_max, sd_u_sum/num_sd_u
endif
if (num_sd_v > 0) then
   write(*,'(9X,A,3ES12.4)') '   V   error SD min/max/mean (m/s): ', sd_v_min, sd_v_max, sd_v_sum/num_sd_v
endif
if (num_sd_r > 0) then
   write(*,'(9X,A,3ES12.4)') 'Radial error SD min/max/mean (m/s): ', sd_r_min, sd_r_max, sd_r_sum/num_sd_r
endif

write(*, '(9X, A)') '+++++++++++++++++'
if (num_sites > 0) then
   site_line = trim(site_names(1))
   do i = 2, num_sites
      site_line = trim(site_line)//', '//trim(site_names(i))
   enddo

   write(*,'(9X,"Radial instruments processed (",I0,"): ",A)') num_sites, trim(site_line)
else
   write(*,'(9X,"Radial instruments processed (0): none")')
endif

end subroutine report_summary


!------------------------------------------------------------
! Clear up memory 
subroutine cleanup()

if (allocated(lon) ) deallocate(lon)
if (allocated(lat) ) deallocate(lat)
if (allocated(u)   ) deallocate(u)
if (allocated(v)   ) deallocate(v)
if (allocated(sd_u)) deallocate(sd_u)
if (allocated(sd_v)) deallocate(sd_v)

if (allocated(rlon)) deallocate(rlon)
if (allocated(rlat)) deallocate(rlat)
if (allocated(rvel)) deallocate(rvel)
if (allocated(rflg)) deallocate(rflg)
if (allocated(rnge)) deallocate(rnge)
if (allocated(rstd)) deallocate(rstd)
if (allocated(phi) ) deallocate(phi)

end subroutine cleanup

end program hf_to_obs

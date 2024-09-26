! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

PROGRAM grid_refl_obs

! Extract reflectivity obs from a DART obs_seq.out file and compute the coordinates of
! these observations on a WRF grid (mass grid points).
!
! David Dowell 26 June 2007
!
! input parameters from command line:
! (1) obs_seq_file  -- path name of DART obs_seq.out file
! (2) refl_min      -- minimum reflectivity threshold for processing observation
! (3) days_begin    -- start of time range to be processed
! (4) seconds_begin -- "                                 "
! (5) days_end      -- end of time range to be processed
! (6) seconds_end   -- "                               "
! (7) wrf_file      -- path name of WRF netcdf file
!
! output:
! (1) text file 

use        types_mod, only : r8, missing_r8, gravity
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_from_key, get_obs_def, get_copy_meta_data, &
                             get_obs_time_range, get_time_range_keys, get_num_obs, &
                             get_next_obs, get_num_times, get_obs_values, init_obs, &
                             get_num_copies, static_init_obs_sequence, &
                             get_qc, destroy_obs_sequence, read_obs_seq_header, & 
                             get_last_obs, destroy_obs, get_num_qc, get_qc_meta_data
use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                             get_obs_def_location,  get_obs_def_type_of_obs
use     obs_kind_mod, only : RADAR_REFLECTIVITY
use        map_utils, only : proj_info, map_init, map_set, latlon_to_ij, &
                             PROJ_LATLON, PROJ_MERC, PROJ_LC, PROJ_PS, &
                             ij_to_latlon, gridwind_to_truewind
use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=),     &
                             is_vertical
use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             set_calendar_type, print_date, GREGORIAN, &
                             operator(*), operator(+), operator(-), &
                             operator(>), operator(<), operator(/), &
                             operator(/=), operator(<=)
use    utilities_mod, only : error_handler, E_ERR, E_MSG, file_exist, &
                             initialize_utilities, finalize_utilities
use    netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! command-line parameters
character(len=129) :: obs_seq_file
real(r8)           :: refl_min
integer            :: seconds_begin
integer            :: days_begin
integer            :: seconds_end
integer            :: days_end
character(len=129) :: wrf_file

! local variables
type(obs_sequence_type)     :: seq
type(obs_def_type)          :: obs_def
type(time_type)             :: beg_time, end_time
type(obs_type)              :: ob
type(location_type)         :: ob_loc

integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id
character(len=129) :: obs_seq_read_format
logical :: pre_I_format
logical :: out_of_range
integer :: key_bounds(2)
integer :: num_obs_in_time_period
integer :: num_refl_obs
integer :: num_refl_obs_in_domain
integer,  allocatable :: keys(:)
integer :: obs_kind_ind
integer :: obs_index
real(r8) :: ob_value(1)

real(r8), allocatable :: lat(:,:)               ! latitude at mass grid points (deg)
real(r8), allocatable :: lon(:,:)               ! longitude at mass grid points (deg)
real(r8), allocatable :: phb(:,:,:)             ! base-state geopotential (m^2 s^-2)
real(r8), allocatable :: ht(:,:,:)              ! height MSL of mass grid points (m)
type(proj_info) :: proj                         ! map projection info.
integer, parameter :: map_sphere = 0, map_lambert = 1, map_polar_stereo = 2, map_mercator = 3
integer map_proj, proj_code
real(r8) :: dx
real(r8) :: stdlon,truelat1,truelat2
real(r8) :: xyz_loc(3)
real(r8) :: iloc, jloc, kloc

real(r8), allocatable :: refl_ob(:,:,:)         ! gridded reflectivity observations (dBZ)

integer               :: bt, sn, we             ! WRF grid dimensions

character(len = 150) :: msgstring
integer               :: i, j, k, o


! netcdf stuff
integer :: var_id, ncid, ierr
character(len=80) :: varname

! command-line parameters stuff
integer :: status, length
character(len=120) :: string


call initialize_utilities('grid_refl_obs')

! Get command-line parameters, using the fortran 2003 intrinsics.

if( COMMAND_ARGUMENT_COUNT() .ne. 7 ) then
  print*, 'INCORRECT # OF ARGUMENTS ON COMMAND LINE:  ', COMMAND_ARGUMENT_COUNT()
  call exit(1)
else

  call GET_COMMAND_ARGUMENT(1,obs_seq_file,length,status)
  if( status .ne. 0 ) then
    print*, 'obs_seq_file NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  endif

  call GET_COMMAND_ARGUMENT(2,string,length,status)
  if( status .ne. 0 ) then
    print*,  'refl_min NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) refl_min
  endif

  call GET_COMMAND_ARGUMENT(3,string,length,status)
  if( status .ne. 0 ) then
    print*,  'days_begin NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) days_begin
  endif

  call GET_COMMAND_ARGUMENT(4,string,length,status)
  if( status .ne. 0 ) then
    print*,  'seconds_begin NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) seconds_begin
  endif

  call GET_COMMAND_ARGUMENT(5,string,length,status)
  if( status .ne. 0 ) then
    print*,  'days_end NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) days_end
  endif

  call GET_COMMAND_ARGUMENT(6,string,length,status)
  if( status .ne. 0 ) then
    print*,  'seconds_end NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  else
    read(string,*) seconds_end
  endif
  
  call GET_COMMAND_ARGUMENT(7,wrf_file,length,status)
  if( status .ne. 0 ) then
    print*, 'wrf_file NOT RETRIEVED FROM COMMAND LINE:  ', status
    call exit(1)
  endif

endif


! Read observations.

call static_init_obs_sequence() ! Initialize the obs sequence module 

obs_seq_file = trim(adjustl(obs_seq_file))
if ( file_exist(obs_seq_file) ) then
  write(msgstring,*)'opening ', obs_seq_file
  call error_handler(E_MSG,'grid_refl_obs',msgstring,source,revision,revdate)
else
  write(msgstring,*)obs_seq_file,&
                    ' does not exist. Finishing up.'
  call error_handler(E_MSG,'grid_refl_obs',msgstring,source,revision,revdate)
  call exit(1)
endif

call read_obs_seq_header(obs_seq_file, &
             num_copies, num_qc, num_obs, max_num_obs, &
             obs_seq_file_id, obs_seq_read_format, pre_I_format, &
             close_the_file = .true.)

print*, 'num_copies = ', num_copies
print*, 'num_qc = ', num_qc
print*, 'num_obs = ', num_obs
print*, 'max_num_obs = ', max_num_obs
print*, 'obs_seq_read_format = ', obs_seq_read_format
print*, 'pre_I_format = ', pre_I_format

call read_obs_seq(obs_seq_file, 0, 0, 0, seq)

MetaDataLoop : do i=1, get_num_copies(seq)
  if(index(get_copy_meta_data(seq,i), 'observation') > 0) obs_index = i
enddo MetaDataLoop


! Open WRF file and obtain grid dimensions.

call check ( nf90_open(wrf_file, NF90_NOWRITE, ncid) )

call check ( nf90_inq_dimid(ncid, 'bottom_top', var_id) )
call check ( nf90_inquire_dimension(ncid, var_id, varname, bt) )

call check ( nf90_inq_dimid(ncid, 'south_north', var_id) )
call check ( nf90_inquire_dimension(ncid, var_id, varname, sn) )

call check ( nf90_inq_dimid(ncid, 'west_east', var_id) )
call check ( nf90_inquire_dimension(ncid, var_id, varname, we) )

call check( nf90_get_att(ncid, nf90_global, 'MAP_PROJ', map_proj) )
call check( nf90_get_att(ncid, nf90_global, 'DX', dx) )
call check( nf90_get_att(ncid, nf90_global, 'TRUELAT1', truelat1) )
call check( nf90_get_att(ncid, nf90_global, 'TRUELAT2', truelat2) )
call check( nf90_get_att(ncid, nf90_global, 'STAND_LON', stdlon) )

! Allocate arrays.

allocate(lat(we,sn))
allocate(lon(we,sn))
allocate(phb(we,sn,bt+1))
allocate(ht(we,sn,bt))
allocate(refl_ob(we,sn,bt))

refl_ob(:,:,:) = missing_r8


! Read WRF grid information.

call check ( nf90_inq_varid(ncid, 'XLAT', var_id))
call check ( nf90_get_var(ncid, var_id, lat, start = (/ 1, 1, 1/)))

call check ( nf90_inq_varid(ncid, 'XLONG', var_id))
call check ( nf90_get_var(ncid, var_id, lon, start = (/ 1, 1, 1/)))

call check ( nf90_inq_varid(ncid, 'PHB', var_id))
call check ( nf90_get_var(ncid, var_id, phb, start = (/ 1, 1, 1, 1/)))

ierr = NF90_close(ncid)


! Set up map projection structure.

call map_init(proj)
if(map_proj == map_sphere) then
  proj_code = PROJ_LATLON
elseif(map_proj == map_lambert) then
  proj_code = PROJ_LC
elseif(map_proj == map_polar_stereo) then
  proj_code = PROJ_PS
elseif(map_proj == map_mercator) then
  proj_code = PROJ_MERC
else
  call error_handler(E_ERR,'grid_refl_obs', &
       'Map projection no supported.', source, revision, revdate)
endif
!call map_set(proj_code,lat(1,1),lon(1,1), &
!             1.0_r8,1.0_r8,dx,stdlon,truelat1,truelat2,proj)
call map_set(proj_code=proj_code, proj=proj, lat1=lat(1,1), lon1=lon(1,1), &
             knowni=1.0_r8, knownj=1.0_r8, dx=dx, stdlon=stdlon, truelat1=truelat1, truelat2=truelat2)

! Compute height (m MSL) of each grid point.

do k=1, bt
  do j=1, sn
    do i=1, we
      ht(i,j,k) = ( phb(i,j,k) + phb(i,j,k+1) ) / (2.0_r8*gravity)
    enddo
  enddo
enddo

! Make sure longitudes are in the range from 0 to 360.

do j=1, sn
  do i=1, we
    do while (lon(i,j) < 0.0_r8)
      lon(i,j) = lon(i,j) + 360.0_r8
    end do
    do while (lon(i,j) > 360.0_r8)
      lon(i,j) = lon(i,j) - 360.0_r8
    end do    
  enddo
enddo


! Process observations, assigning to grid points any reflectivity observations that lie within
! the domain and specified time range.

beg_time = set_time(seconds_begin, days_begin)
end_time = set_time(seconds_end, days_end)
call get_obs_time_range(seq, beg_time, end_time, key_bounds, &
                        num_obs_in_time_period, out_of_range)

print*, 'total number of observations in time period = ', num_obs_in_time_period

allocate(keys(num_obs_in_time_period))
call get_time_range_keys(seq, key_bounds, num_obs_in_time_period, keys)

num_refl_obs = 0
num_refl_obs_in_domain = 0

do o = 1, num_obs_in_time_period

  call get_obs_from_key(seq, keys(o), ob)
  call get_obs_def(ob, obs_def)
  ob_loc = get_obs_def_location(obs_def)
  obs_kind_ind = get_obs_def_type_of_obs(obs_def)

  if ( (obs_kind_ind == RADAR_REFLECTIVITY) .and. (is_vertical(ob_loc, "HEIGHT")) ) then
  
    num_refl_obs = num_refl_obs + 1

!   xyz_loc(1) is longitude,  xyz_loc(2) is latitude,  xyz_loc(3) is height
    xyz_loc = get_location(ob_loc)
    call latlon_to_ij(proj,xyz_loc(2),xyz_loc(1),iloc,jloc)
    if ( (iloc >= 1 .and. iloc <= we .and. jloc >= 1 .and. jloc <= sn) ) then

      i = nint(iloc)
      j = nint(jloc)
      
      kloc = 0.0
      do k=1, bt-1
        if ( (xyz_loc(3).ge.ht(i,j,k)) .and. (xyz_loc(3).le.ht(i,j,k+1)) ) then
          kloc = real(k) + (xyz_loc(3)-ht(i,j,k)) / (ht(i,j,k+1)-ht(i,j,k))
        endif
      enddo
      k = nint(kloc)

      if (k.ne.0) then

        num_refl_obs_in_domain = num_refl_obs_in_domain + 1
        call get_obs_values(ob, ob_value, obs_index)
        if ( (refl_ob(i,j,k).eq.missing_r8) .or. (ob_value(1).gt.refl_ob(i,j,k)) ) then
          refl_ob(i,j,k) = ob_value(1)
        endif

      endif

    endif

  end if       ! if (obs_kind_ind == RADAR_REFLECTIVITY)

enddo        ! do o = 1, num_obs_in_time_period

print*, 'number of reflectivity observations in time period = ', num_refl_obs
print*, 'number of reflectivity observations in time period and domain = ', num_refl_obs_in_domain


! Output gridded reflectivity values that exceed threshold.

open(unit=11, file='refl_obs.txt', status='unknown')
do k=1, bt
  do j=1, sn
    do i=1, we
      if ( (refl_ob(i,j,k).ne.missing_r8) .and. (refl_ob(i,j,k).gt.refl_min) ) then
        write(11,*) i, j, k, refl_ob(i,j,k)
      endif
    enddo
  enddo
enddo
close(11)


! Deallocate arrays.

deallocate(lat)
deallocate(lon)
deallocate(phb)
deallocate(ht)
deallocate(refl_ob)
deallocate(keys)

call finalize_utilities('grid_refl_obs')

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'grid_refl_obs', &
         trim(nf90_strerror(istatus)), source, revision, revdate)
  end subroutine check

END PROGRAM grid_refl_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program wrf_dart_obs_preprocess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   wrf_dart_obs_preprocess - WRF-DART utility program that at a
!                             minimum removes all observations outside
!                             of a WRF domain and add observations from
!                             supplimental obs sequences.  The program
!                             assumes all data is from one observation
!                             time.  In addition, this program allows 
!                             the user to do the following functions:
!
!     - remove observations near the lateral boundaries
!     - increase observation error near lateral boundaries
!     - remove observations above certain pressure/height levels
!     - remove observations where the model and obs topography are large
!     - remove significant level rawinsonde data
!     - remove rawinsonde observations near TC core
!     - superob aircraft and satellite wind data
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use        types_mod, only : r8, i8
use obs_sequence_mod, only : obs_sequence_type, static_init_obs_sequence, &
                             read_obs_seq_header, destroy_obs_sequence, &
                             get_num_obs, write_obs_seq 
use    utilities_mod, only : find_namelist_in_file, check_namelist_read, &
                             initialize_utilities, finalize_utilities
use  netcdf_utilities_mod, only : nc_check
use     obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT, ACARS_U_WIND_COMPONENT, &
                             MARINE_SFC_U_WIND_COMPONENT, LAND_SFC_U_WIND_COMPONENT, &
                             METAR_U_10_METER_WIND, GPSRO_REFRACTIVITY, &
                             SAT_U_WIND_COMPONENT, PROFILER_U_WIND_COMPONENT, VORTEX_LAT
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time
use ensemble_manager_mod, only : ensemble_type, init_ensemble_manager, end_ensemble_manager
use        model_mod, only : static_init_model
use           netcdf

implicit none

! ----------------------------------------------------------------------
! Declare namelist parameters
! ----------------------------------------------------------------------

!  Generic parameters
character(len=129) :: file_name_input    = 'obs_seq.old',        &
                      file_name_output   = 'obs_seq.new',        &
                      sonde_extra        = 'obs_seq.rawin',      &
                      acars_extra        = 'obs_seq.acars',      &
                      land_sfc_extra     = 'obs_seq.land_sfc',   &
                      metar_extra        = 'obs_seq.metar',      &
                      marine_sfc_extra   = 'obs_seq.marine',     &
                      sat_wind_extra     = 'obs_seq.satwnd',     &
                      profiler_extra     = 'obs_seq.profiler',   &
                      gpsro_extra        = 'obs_seq.gpsro',      &
                      trop_cyclone_extra = 'obs_seq.tc'
integer            :: max_num_obs              = 100000   ! Largest number of obs in one sequence
logical            :: overwrite_obs_time       = .false.  ! true to overwrite all observation times

!  boundary-specific parameters
real(r8)           :: obs_boundary             = 0.0_r8   ! number of grid points to remove obs near boundary
logical            :: increase_bdy_error       = .false.  ! true to increase obs error near boundary
real(r8)           :: maxobsfac                = 2.5_r8   ! maximum increase in obs error near boundary
real(r8)           :: obsdistbdy               = 15.0_r8  ! number of grid points to increase obs err.

!  parameters used to reduce observations
logical            :: sfc_elevation_check      = .false.  ! remove obs where model-obs topography is large
real(r8)           :: sfc_elevation_tol        = 300.0_r8  ! largest difference between model and obs. topo.
real(r8)           :: obs_pressure_top         = 0.0_r8  ! remove all obs at lower pressure
real(r8)           :: obs_height_top           = 2.0e10_r8  ! remove all obs at higher height

!  Rawinsonde-specific parameters
logical            :: include_sig_data         = .true.   ! include significant-level data
real(r8)           :: tc_sonde_radii           = -1.0_r8  ! remove sonde obs closer than this to TC

!  aircraft-specific parameters
logical            :: superob_aircraft         = .false.  ! super-ob aircraft data
real(r8)           :: aircraft_horiz_int       = 36.0_r8  ! horizontal interval for super-ob (km)
real(r8)           :: aircraft_pres_int        = 2500.0_r8  ! pressure interval for super-ob

!  sat wind specific parameters
logical            :: superob_sat_winds        = .false.    ! super-ob sat wind data
real(r8)           :: sat_wind_horiz_int       = 100.0_r8   ! horizontal interval for super-ob
real(r8)           :: sat_wind_pres_int        = 2500.0_r8  ! pressure interval for super-ob
logical            :: overwrite_ncep_satwnd_qc = .false.    ! true to overwrite NCEP QC (see instructions)

!  surface obs. specific parameters
logical            :: overwrite_ncep_sfc_qc    = .false.  ! true to overwrite NCEP QC (see instructions)

namelist /wrf_obs_preproc_nml/file_name_input, file_name_output,      &
         include_sig_data, superob_aircraft, superob_sat_winds,     &
         sfc_elevation_check, overwrite_ncep_sfc_qc, overwrite_ncep_satwnd_qc, &
         aircraft_pres_int, sat_wind_pres_int, sfc_elevation_tol,   & 
         obs_pressure_top, obs_height_top, obs_boundary, sonde_extra, metar_extra,   &
         acars_extra, land_sfc_extra, marine_sfc_extra, sat_wind_extra, profiler_extra, &
         trop_cyclone_extra, gpsro_extra, tc_sonde_radii, increase_bdy_error,      &
         maxobsfac, obsdistbdy, sat_wind_horiz_int, aircraft_horiz_int, &
         overwrite_obs_time

! ----------------------------------------------------------------------
! Declare other variables
! ----------------------------------------------------------------------

character(len=129)      :: obs_seq_read_format
character(len=80)       :: name

integer                 :: io, iunit, fid, var_id, obs_seq_file_id, num_copies, &
                           num_qc, num_obs, max_obs_seq, nx, ny, gday, gsec

real(r8)                :: real_nx, real_ny

logical                 :: file_exist, pre_I_format

type(obs_sequence_type) :: seq_all, seq_rawin, seq_sfc, seq_acars, seq_satwnd, &
                           seq_prof, seq_tc, seq_gpsro, seq_other

type(time_type)         :: anal_time

type(ensemble_type)     :: dummy_ens

call initialize_utilities("wrf_dart_obs_preprocess")

print*,'Enter target assimilation time (gregorian day, second): '
read*,gday,gsec
call set_calendar_type(GREGORIAN)
anal_time = set_time(gsec, gday)

call static_init_obs_sequence()
call static_init_model()
call init_ensemble_manager(dummy_ens, 1, 1_i8)

call find_namelist_in_file("input.nml", "wrf_obs_preproc_nml", iunit)
read(iunit, nml = wrf_obs_preproc_nml, iostat = io)
call check_namelist_read(iunit, io, "wrf_obs_preproc_nml")

!  open a wrfinput file, which is on this domain
call nc_check( nf90_open(path = "wrfinput_d01", mode = nf90_nowrite, ncid = fid), &
               'main', 'open wrfinput_d01' )
call nc_check( nf90_inq_dimid(fid, "west_east", var_id), &
               'main', 'inq. dimid west_east' )
call nc_check( nf90_inquire_dimension(fid, var_id, name, nx), &
               'main', 'inquire dimension west_east' )
call nc_check( nf90_inq_dimid(fid, "south_north", var_id), &
               'main', 'inq. dimid south_north' )
call nc_check( nf90_inquire_dimension(fid, var_id, name, ny), &
               'main', 'inquire dimension south_north' )
call nc_check( nf90_close(fid), 'main', 'close wrfinput_d01' )

! several places need a real(r8) version of nx and ny, so set them up
! here so they're ready to use.  previous versions of this code used
! the conversions dble(nx) but that doesn't work if you are compiling
! this code with r8 redefined to be r4.  on some platforms it corrupts
! the arguments to the function call.  these vars are guarenteed to be
! the right size for code compiled with reals as either r8 or r4.
real_nx = nx
real_ny = ny

!  if obs_seq file exists, read in the data, otherwise, create a blank one.
inquire(file = trim(adjustl(file_name_input)), exist = file_exist)
if ( file_exist ) then

  call read_obs_seq_header(file_name_input, num_copies, num_qc, num_obs, max_obs_seq, &
        obs_seq_file_id, obs_seq_read_format, pre_I_format, close_the_file = .true.)

else

  num_copies = 1  ;  num_qc = 1  ;  max_obs_seq = max_num_obs * 3
  call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_all)

end if

!  create obs sequences for different obs types
call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_rawin)
call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_sfc)
call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_acars)
call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_satwnd)
call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_prof)
call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_gpsro)
call create_new_obs_seq(num_copies, num_qc, 100,         seq_tc)
call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_other)

!  read input obs_seq file, divide into platforms
call read_and_parse_input_seq(file_name_input, real_nx, real_ny, obs_boundary, &
include_sig_data, obs_pressure_top, obs_height_top, sfc_elevation_check, &
sfc_elevation_tol, overwrite_ncep_sfc_qc, overwrite_ncep_satwnd_qc, &
overwrite_obs_time, anal_time, seq_rawin, seq_sfc, seq_acars, seq_satwnd, & 
seq_tc, seq_gpsro, seq_other)

!  add supplimental rawinsonde observations from file
call add_supplimental_obs(sonde_extra, seq_rawin, max_obs_seq, &
RADIOSONDE_U_WIND_COMPONENT, nx, ny, obs_boundary, include_sig_data, &
obs_pressure_top, obs_height_top, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time)

!  add supplimental ACARS observations from file
call add_supplimental_obs(acars_extra, seq_acars, max_obs_seq, &
ACARS_U_WIND_COMPONENT, nx, ny, obs_boundary, include_sig_data, &
obs_pressure_top, obs_height_top, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time)

!  add supplimental marine observations from file
call add_supplimental_obs(marine_sfc_extra, seq_sfc, max_obs_seq, &
MARINE_SFC_U_WIND_COMPONENT, nx, ny, obs_boundary, include_sig_data, &
obs_pressure_top, obs_height_top, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time)

!  add supplimental land surface observations from file
call add_supplimental_obs(land_sfc_extra, seq_sfc, max_obs_seq, &
LAND_SFC_U_WIND_COMPONENT, nx, ny, obs_boundary, include_sig_data, &
obs_pressure_top, obs_height_top, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time)

!  add supplimental metar observations from file
call add_supplimental_obs(metar_extra, seq_sfc, max_obs_seq, &
METAR_U_10_METER_WIND, nx, ny, obs_boundary, include_sig_data, &
obs_pressure_top, obs_height_top, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time)

!  add supplimental satellite wind observations from file
call add_supplimental_obs(sat_wind_extra, seq_satwnd, max_obs_seq, &
SAT_U_WIND_COMPONENT, nx, ny, obs_boundary, include_sig_data, &
obs_pressure_top, obs_height_top, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time)

!  add supplimental profiler observations from file
call add_supplimental_obs(profiler_extra, seq_prof, max_obs_seq, &
PROFILER_U_WIND_COMPONENT, nx, ny, obs_boundary, include_sig_data, &
obs_pressure_top, obs_height_top, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time)

!  add supplimental GPSRO observations from file
call add_supplimental_obs(gpsro_extra, seq_gpsro, max_obs_seq, &
GPSRO_REFRACTIVITY, nx, ny, obs_boundary, include_sig_data, &
obs_pressure_top, obs_height_top, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time)

!  add supplimental tropical cyclone vortex observations from file
call add_supplimental_obs(trop_cyclone_extra, seq_tc, max_obs_seq, &
VORTEX_LAT, nx, ny, obs_boundary, include_sig_data, &
obs_pressure_top, obs_height_top, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time)

!  remove all sonde observations within radius of TC if desired
if ( tc_sonde_radii > 0.0_r8 ) call remove_sondes_near_tc(seq_tc, & 
                                               seq_rawin, tc_sonde_radii)

!  super-ob ACARS data
if ( superob_aircraft ) call superob_aircraft_data(seq_acars, anal_time, &
                                     aircraft_horiz_int, aircraft_pres_int)

!  super-ob satellite wind data
if ( superob_sat_winds ) call superob_sat_wind_data(seq_satwnd, anal_time, &
                                     sat_wind_horiz_int, sat_wind_pres_int)

max_obs_seq = get_num_obs(seq_tc)     + get_num_obs(seq_rawin) + &
              get_num_obs(seq_sfc)    + get_num_obs(seq_acars) + &
              get_num_obs(seq_satwnd) + get_num_obs(seq_prof)  + &
              get_num_obs(seq_gpsro)  +  get_num_obs(seq_other)

call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_all)

call build_master_sequence(seq_tc, seq_all)
call destroy_obs_sequence(seq_tc)

call build_master_sequence(seq_rawin, seq_all)
call destroy_obs_sequence(seq_rawin)

call build_master_sequence(seq_sfc, seq_all)
call destroy_obs_sequence(seq_sfc)

call build_master_sequence(seq_acars, seq_all)
call destroy_obs_sequence(seq_acars)

call build_master_sequence(seq_gpsro, seq_all)
call destroy_obs_sequence(seq_gpsro)

call build_master_sequence(seq_satwnd, seq_all)
call destroy_obs_sequence(seq_satwnd)

call build_master_sequence(seq_prof, seq_all)
call destroy_obs_sequence(seq_prof)

call build_master_sequence(seq_other, seq_all)
call destroy_obs_sequence(seq_other)

write(6,*) 'Total number of observations:', get_num_obs(seq_all)

!  increase the observation error along the lateral boundary
if ( increase_bdy_error ) call increase_obs_err_bdy(seq_all, &
                              obsdistbdy, maxobsfac, real_nx, real_ny)

!  write the observation sequence to file
call write_obs_seq(seq_all, file_name_output)
call destroy_obs_sequence(seq_all)

call finalize_utilities("wrf_dart_obs_preprocess")

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   aircraft_obs_check - function that determines whether to include an
!                     aircraft observation in the sequence.  For now,
!                     this function is a placeholder and returns true.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function aircraft_obs_check()

use     types_mod, only : r8

implicit none

logical  :: aircraft_obs_check

aircraft_obs_check = .true.

return
end function aircraft_obs_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   add_supplimental_obs - subroutine that reads observation data from
!                          a supplimental obs sequence file, performs
!                          validation checks and adds it to the
!                          platform-specific obs sequence.
!
!    filename    - name of supplimental obs sequence file
!    obs_seq     - platform-specific obs sequence
!    max_obs_seq - maximum number of observations in sequence
!    plat_kind   - integer kind of platform (used for print statements)
!    nx          - number of grid points in x direction
!    ny          - number of grid points in y direction
!    obs_bdy     - grid point buffer to remove observations
!    siglevel    - true to include sonde significant level data
!    ptop        - lowest pressure to include in sequence
!    htop        - highest height level to include in sequence
!    sfcelev     - true to perform surface obs. elevation check
!    elev_max    - maximum difference between model and obs. height
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_supplimental_obs(filename, obs_seq, max_obs_seq, plat_kind, &
                                 nx, ny, obs_bdy, siglevel, ptop, htop, &
                                 sfcelev, elev_max, overwrite_time, atime)

use         types_mod, only : r8
use  time_manager_mod, only : time_type, operator(>=)
use      location_mod, only : location_type, get_location, is_vertical
use  obs_sequence_mod, only : obs_sequence_type, obs_type, init_obs, set_obs_def, &
                              get_num_copies, get_num_qc, read_obs_seq, copy_obs, &
                              get_first_obs, get_obs_def, get_next_obs, &
                              get_last_obs, insert_obs_in_seq, destroy_obs_sequence
use       obs_def_mod, only : obs_def_type, get_obs_def_type_of_obs, set_obs_def_time, &
                              get_obs_def_location, get_obs_def_time
use      obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT, ACARS_U_WIND_COMPONENT, &
                              LAND_SFC_U_WIND_COMPONENT, MARINE_SFC_U_WIND_COMPONENT, &
                              METAR_U_10_METER_WIND, GPSRO_REFRACTIVITY, &
                              SAT_U_WIND_COMPONENT, VORTEX_LAT
use         model_mod, only : get_domain_info 

implicit none

character(len=129), intent(in)         :: filename
type(time_type),         intent(in)    :: atime
type(obs_sequence_type), intent(inout) :: obs_seq
integer, intent(in)                    :: max_obs_seq, plat_kind, nx, ny
logical, intent(in)                    :: siglevel, sfcelev, overwrite_time
real(r8), intent(in)                   :: obs_bdy, ptop, htop, elev_max

integer  :: nloc, okind, dom_id
logical  :: file_exist, last_obs, pass_checks, first_obs
real(r8) :: xyz_loc(3), xloc, yloc
real(r8) :: real_nx, real_ny


type(location_type)     :: obs_loc_list(max_obs_seq), obs_loc
type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: supp_obs_seq
type(obs_type)          :: obs_in, prev_obsi, prev_obso, obs
type(time_type)         :: obs_time, prev_time

inquire(file = trim(adjustl(filename)), exist = file_exist)
if ( .not. file_exist )  return

! see comment in main routine about why these are needed
real_nx = nx
real_ny = ny

select case (plat_kind)

  case (RADIOSONDE_U_WIND_COMPONENT)
    write(6,*) 'Adding Supplimental Rawinsonde Data'
  case (ACARS_U_WIND_COMPONENT)
    write(6,*) 'Adding Supplimental ACARS Data'
  case (MARINE_SFC_U_WIND_COMPONENT)
    write(6,*) 'Adding Supplimental Marine Surface Data'
  case (LAND_SFC_U_WIND_COMPONENT)
    write(6,*) 'Adding Supplimental Land Surface Data'
  case (METAR_U_10_METER_WIND)
    write(6,*) 'Adding Supplimental METAR Data'
  case (SAT_U_WIND_COMPONENT)
    write(6,*) 'Adding Supplimental Satellite Wind Data'
  case (VORTEX_LAT)
    write(6,*) 'Adding Supplimental Tropical Cyclone Data'
  case (GPSRO_REFRACTIVITY)
    write(6,*) 'Adding Supplimental GPS RO Data'

end select

call init_obs(obs_in,    get_num_copies(obs_seq), get_num_qc(obs_seq))
call init_obs(obs,       get_num_copies(obs_seq), get_num_qc(obs_seq))
call init_obs(prev_obsi, get_num_copies(obs_seq), get_num_qc(obs_seq))
call init_obs(prev_obso, get_num_copies(obs_seq), get_num_qc(obs_seq))

!  create list of observations in plaform sequence
call build_obs_loc_list(obs_seq, max_obs_seq, nloc, obs_loc_list)

!  find the last observation in the sequence
if ( get_last_obs(obs_seq, prev_obso) ) then

  first_obs = .false.
  call get_obs_def(prev_obso, obs_def)
  prev_time = get_obs_def_time(obs_def)

else

  first_obs = .true.

end if

last_obs = .false.
call read_obs_seq(trim(adjustl(filename)), 0, 0, 0, supp_obs_seq)
if ( .not. get_first_obs(supp_obs_seq, obs_in) ) last_obs = .true.

ObsLoop:  do while ( .not. last_obs ) ! loop over all observations in a sequence

  !  read data from observation
  call get_obs_def(obs_in, obs_def)
  okind   = get_obs_def_type_of_obs(obs_def)
  obs_loc = get_obs_def_location(obs_def)
  xyz_loc = get_location(obs_loc)
  call get_domain_info(xyz_loc(1),xyz_loc(2),dom_id,xloc,yloc)

  !  check if the observation is within the domain
  if ( ((xloc < (obs_bdy+1.0_r8) .or. xloc > (real_nx-obs_bdy-1.0_r8) .or. &
         yloc < (obs_bdy+1.0_r8) .or. yloc > (real_ny-obs_bdy-1.0_r8)) .and. &
         (dom_id == 1)) .or. dom_id < 1 ) then

    prev_obsi = obs_in
    call get_next_obs(supp_obs_seq, prev_obsi, obs_in, last_obs)
    cycle ObsLoop

  end if

  !  check if the observation is within vertical bounds of domain
  if ( (is_vertical(obs_loc, "PRESSURE") .and. xyz_loc(3) < ptop) .or. &
       (is_vertical(obs_loc, "HEIGHT")   .and. xyz_loc(3) > htop) ) then

    prev_obsi = obs_in
    call get_next_obs(supp_obs_seq, prev_obsi, obs_in, last_obs)
    cycle ObsLoop

  end if

  !  check if the observation already exists
  if ( .not. original_observation(obs_loc, obs_loc_list, nloc) ) then

    prev_obsi = obs_in
    call get_next_obs(supp_obs_seq, prev_obsi, obs_in, last_obs)
    cycle ObsLoop

  end if

  !  overwrite the observation time with the analysis time if desired
  if ( overwrite_time ) then

    call set_obs_def_time(obs_def, atime)
    call set_obs_def(obs_in, obs_def)

  end if

  ! perform platform-specific checks
  select case (plat_kind)

    case (RADIOSONDE_U_WIND_COMPONENT)
      pass_checks = rawinsonde_obs_check(obs_loc, okind, siglevel, &
                                                sfcelev, elev_max)
    case (ACARS_U_WIND_COMPONENT)
      pass_checks = aircraft_obs_check()
    case (MARINE_SFC_U_WIND_COMPONENT)
      pass_checks = surface_obs_check(sfcelev, elev_max, xyz_loc)
    case (LAND_SFC_U_WIND_COMPONENT)
      pass_checks = surface_obs_check(sfcelev, elev_max, xyz_loc)
    case (METAR_U_10_METER_WIND)
      pass_checks = surface_obs_check(sfcelev, elev_max, xyz_loc)
    case (SAT_U_WIND_COMPONENT)
      pass_checks = sat_wind_obs_check()
    case default
      pass_checks = .true.

  end select

  if ( pass_checks ) then

    call copy_obs(obs, obs_in)
    call get_obs_def(obs, obs_def)
    obs_time = get_obs_def_time(obs_def)

    if (obs_time >= prev_time .and. (.not. first_obs)) then  ! same time or later than previous obs
      call insert_obs_in_seq(obs_seq, obs, prev_obso)
    else                                                     ! earlier, search from start of seq
      call insert_obs_in_seq(obs_seq, obs)
    end if

    first_obs = .false.
    prev_obso = obs
    prev_time = obs_time

  end if

  prev_obsi = obs_in
  call get_next_obs(supp_obs_seq, prev_obsi, obs_in, last_obs)

end do ObsLoop

call destroy_obs_sequence(supp_obs_seq)

return
end subroutine add_supplimental_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_new_obs_seq - subroutine that is used to create a new 
!                        observation sequence.
!
!    num_copies - number of copies associated with each observation
!    num_qc     - number of quality control reports in each obs.
!    max_num    - maximum number of observations in sequence
!    seq        - observation sequence
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_new_obs_seq(num_copies, num_qc, max_num, seq)

use obs_sequence_mod, only : obs_sequence_type, init_obs_sequence, &
                             set_copy_meta_data, set_qc_meta_data

implicit none

integer, intent(in) :: num_copies, num_qc, max_num
type(obs_sequence_type), intent(out) :: seq

character(len=129) :: copy_meta_data, qc_meta_data
integer :: i

call init_obs_sequence(seq, num_copies, num_qc, max_num)
do i = 1, num_copies
   copy_meta_data = 'NCEP BUFR observation'
   call set_copy_meta_data(seq, i, copy_meta_data)
end do
do i = 1, num_qc
   qc_meta_data = 'NCEP QC index'
   call set_qc_meta_data(seq, i, qc_meta_data)
end do

return
end subroutine create_new_obs_seq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   build_master_sequence - subroutine used to take observations from
!                           a smaller observation sequence and appends
!                           them to a larger observation sequence.
!                           Note that this routine only works if the
!                           observations are at the same time.
!
!    seq_type - observation sequence with one observation type
!    seq_all  - observation sequence with more observations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_master_sequence(seq_type, seq_all)

use time_manager_mod, only : time_type, operator(>=)
use obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs, & 
                             get_first_obs, copy_obs, insert_obs_in_seq, & 
                             get_next_obs, get_obs_def, get_num_copies, &
                             get_num_qc, get_obs_def
use      obs_def_mod, only : obs_def_type, get_obs_def_time

implicit none

type(obs_sequence_type), intent(in)    :: seq_type
type(obs_sequence_type), intent(inout) :: seq_all

logical :: last_obs, first_obs
type(obs_def_type) :: obs_def
type(obs_type)     :: obs_in, obs, prev_obsi, prev_obsa
type(time_type)    :: obs_time, prev_time

last_obs = .false.  ;  first_obs = .true.
call init_obs(obs_in,    get_num_copies(seq_type), get_num_qc(seq_type))
call init_obs(obs,       get_num_copies(seq_type), get_num_qc(seq_type))
call init_obs(prev_obsi, get_num_copies(seq_type), get_num_qc(seq_type))
call init_obs(prev_obsa, get_num_copies(seq_type), get_num_qc(seq_type))

if ( .not. get_first_obs(seq_type, obs_in) )  return

do while ( .not. last_obs )

  call copy_obs(obs, obs_in)
  call get_obs_def(obs, obs_def)
  obs_time = get_obs_def_time(obs_def)

  if (obs_time >= prev_time .and. (.not. first_obs)) then  ! same time or later than previous obs
    call insert_obs_in_seq(seq_all, obs, prev_obsa) 
  else                                                      ! earlier, search from start of seq
    call insert_obs_in_seq(seq_all, obs)
  end if

  first_obs = .false.
  prev_obsi = obs_in
  prev_obsa = obs
  prev_time = obs_time

  call get_next_obs(seq_type, prev_obsi, obs_in, last_obs)

end do

return
end subroutine build_master_sequence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   build_obs_loc_list - subroutine that creates an array of locations
!                        of the observations in a sequence.
!
!    obs_seq      - observation sequence to read locations from
!    maxobs       - maximum number of observations in a sequence
!    nloc         - number of individual locations
!    obs_loc_list - array of observation locations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine build_obs_loc_list(seq, maxobs, nloc, obs_loc_list)

use       location_mod, only : location_type
use   obs_sequence_mod, only : obs_sequence_type, get_num_copies, &
                               get_num_qc, init_obs, get_first_obs, &
                               get_obs_def, get_next_obs, obs_type
use        obs_def_mod, only : obs_def_type, get_obs_def_location

implicit none

integer, intent(in)                 :: maxobs
type(obs_sequence_type), intent(in) :: seq
integer, intent(out)                :: nloc 
type(location_type), intent(out)    :: obs_loc_list(maxobs)

logical             :: last_obs
type(obs_type)      :: obs, prev_obs
type(obs_def_type)  :: obs_def
type(location_type) :: obs_loc

call init_obs(obs,      get_num_copies(seq), get_num_qc(seq))
call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))

last_obs = .false.  ;  nloc = 0
if ( .not. get_first_obs(seq, obs) ) last_obs = .true.

do while ( .not. last_obs ) ! loop over all observations in a sequence

  call get_obs_def(obs, obs_def)
  obs_loc = get_obs_def_location(obs_def)

  !  construct a list of observation locations
  if ( original_observation(obs_loc, obs_loc_list, nloc) ) then

    nloc = nloc + 1
    obs_loc_list(nloc) = obs_loc

  end if
  prev_obs = obs
  call get_next_obs(seq, prev_obs, obs, last_obs)

end do

return
end subroutine build_obs_loc_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_obs_type - subroutine that is used to create an observation 
!                     type from observation data.
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vloc  - vertical location of observation
!    vcord - DART vertical coordinate integer
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error
!    day   - gregorian day of the observation
!    sec   - gregorian second of the observation
!    qc    - integer quality control value
!    obs   - observation type that includes the observation information
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_obs_type(lat, lon, vloc, vcord, obsv, okind, oerr, qc, otime, obs)

use types_mod,        only : r8
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location
use     location_mod, only : location_type, set_location
use time_manager_mod, only : time_type

implicit none

integer, intent(in)           :: okind, vcord
real(r8), intent(in)          :: lat, lon, vloc, obsv, oerr, qc
type(time_type), intent(in)   :: otime
type(obs_type), intent(inout) :: obs

real(r8)              :: obs_val(1), qc_val(1)
type(obs_def_type)    :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vloc, vcord))
call set_obs_def_type_of_obs(obs_def, okind)
call set_obs_def_time(obs_def, otime)
call set_obs_def_error_variance(obs_def, oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

return
end subroutine create_obs_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   increase_obs_err_bdy - subroutine that increases the observation
!                          error based on proximity to the lateral
!                          boundary.
!
!    seq    - observation sequence
!    obsbdy - number of grid points near boundary to increase error
!    maxfac - factor to increase observation error at boundary
!    nx     - number of grid points in the x direction
!    ny     - number of grid points in the y direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine increase_obs_err_bdy(seq, obsbdy, maxfac, nx, ny)

use         types_mod, only : r8
use      location_mod, only : location_type, get_location
use  obs_sequence_mod, only : obs_sequence_type, obs_type, init_obs, &
                              get_num_copies, get_num_qc, get_first_obs, &
                              get_obs_def, set_obs_def, set_obs, &
                              get_next_obs, get_obs_key
use       obs_def_mod, only : obs_def_type, get_obs_def_error_variance, &
                              set_obs_def_error_variance, get_obs_def_location
use         model_mod, only : get_domain_info

implicit none

type(obs_sequence_type), intent(inout) :: seq
real(r8), intent(in)                   :: obsbdy, maxfac, nx, ny

integer            :: dom_id
logical            :: last_obs
real(r8)           :: mobse, bobse, xyz_loc(3), xloc, yloc, bdydist, obsfac

type(obs_def_type) :: obs_def
type(obs_type)     :: obs, prev_obs

write(6,*) 'Increasing the Observation Error Near the Lateral Boundary'

call init_obs(obs,      get_num_copies(seq), get_num_qc(seq))
call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))

! compute slope and intercept for error increase factor
mobse = (maxfac - 1.0_r8) / (1.0_r8 - obsbdy)
bobse = maxfac - mobse

last_obs = .false.
if ( .not. get_first_obs(seq, obs) ) last_obs = .true.

do while ( .not. last_obs )

  !  get location information relative to domain 1 (skip nests)
  call get_obs_def(obs, obs_def)
  xyz_loc = get_location(get_obs_def_location(obs_def))
  call get_domain_info(xyz_loc(1),xyz_loc(2),dom_id,xloc,yloc,1)

  !  compute distance to boundary, increase based on this distance
  bdydist = min(xloc-1.0_r8, yloc-1.0_r8, nx-xloc, ny-yloc)
  if ( bdydist <= obsbdy ) then

    obsfac = mobse * bdydist + bobse
    call set_obs_def_error_variance(obs_def, &
           get_obs_def_error_variance(obs_def) * obsfac * obsfac)
    call set_obs_def(obs, obs_def)
    call set_obs(seq, obs, get_obs_key(obs))

  end if
  prev_obs = obs
  call get_next_obs(seq, prev_obs, obs, last_obs)

end do

return
end subroutine increase_obs_err_bdy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   isManLevel - function that returns a logical true if the input 
!                pressure level is a mandatory rawinsonde level.
!
!    plevel - pressure level to check (Pa)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function isManLevel(plevel)

use types_mod,        only : r8

implicit none

real(r8), intent(in) :: plevel

integer, parameter :: nman = 16
integer :: kk
logical :: isManLevel
real(r8) :: raw_man_levels(nman) =  (/  &
       100000.0_r8, 92500.0_r8, 85000.0_r8, 70000.0_r8, 50000.0_r8, 40000.0_r8, &
        30000.0_r8, 25000.0_r8, 20000.0_r8, 15000.0_r8, 10000.0_r8,  7000.0_r8, &
         5000.0_r8,  3000.0_r8,  2000.0_r8,  1000.0_r8 /)

isManLevel = .false.
do kk = 1, nman
  if ( plevel == raw_man_levels(kk) ) then
    isManLevel = .true.
    return 
  end if
end do

return
end function isManLevel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   original_observation - function that returns true if the location
!                          is not within an array of locations
!
!    obsloc      - location to check
!    obsloc_list - array of locations to look through
!    nloc        - number of locations in array
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function original_observation(obsloc, obsloc_list, nloc)

use     types_mod, only : r8
use  location_mod, only : location_type, get_dist

real(r8), parameter :: dist_epsilon = 0.00001_r8

integer, intent(in)             :: nloc
type(location_type), intent(in) :: obsloc, obsloc_list(nloc)

integer :: n
logical :: original_observation

original_observation = .true.

do n = 1, nloc

  if ( get_dist(obsloc, obsloc_list(n), 1, 1, .true.) <= dist_epsilon ) then
    original_observation = .false.
    return
  end if

end do

return
end function original_observation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   rawinsonde_obs_check - function that performs obsrvation checks
!                          specific to rawinsonde observations.
!
!    obs_loc    - observation location
!    obs_kind   - DART observation kind
!    siglevel   - true to include significant level data
!    elev_check - true to check differene between model and obs elev.
!    elev_max   - maximum difference between model and obs elevation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rawinsonde_obs_check(obs_loc, obs_kind, siglevel, &
                                elev_check, elev_max)

use     types_mod, only : r8
use  obs_kind_mod, only : RADIOSONDE_SURFACE_ALTIMETER, QTY_SURFACE_ELEVATION
use     model_mod, only : model_interpolate
use  location_mod, only : location_type, set_location, get_location

implicit none

type(location_type), intent(in) :: obs_loc
integer, intent(in)             :: obs_kind
logical, intent(in)             :: siglevel, elev_check
real(r8), intent(in)            :: elev_max

integer  :: istatus(1)
logical  :: rawinsonde_obs_check
real(r8) :: xyz_loc(3), xmod(1), hsfc(1)

rawinsonde_obs_check = .true.
xyz_loc = get_location(obs_loc)

if ( obs_kind /= RADIOSONDE_SURFACE_ALTIMETER ) then

  !  check if vertical level is mandatory level
  if ( (.not. siglevel) .and. (.not. isManLevel(xyz_loc(3))) ) then
    rawinsonde_obs_check = .false.
    return
  end if

else

  !  perform elevation check for altimeter
  if ( elev_check ) then

    call model_interpolate(dummy_ens, 1, obs_loc, QTY_SURFACE_ELEVATION, hsfc, istatus)
    if ( abs(hsfc(1) - xyz_loc(3)) > elev_max ) rawinsonde_obs_check = .false.

  end if

end if

return
end function rawinsonde_obs_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   read_and_parse_input_seq - subroutine that reads a generic
!                              observation sequence and divides the
!                              obs into sequences for each platform.
!
!    filename      - name of input obs sequence
!    nx            - number of grid points in x direction
!    ny            - number of grid points in y direction
!    obs_bdy       - grid point buffer to remove observations
!    siglevel      - true to include sonde significant level data
!    ptop          - lowest pressure to include in sequence
!    htop          - highest height level to include in sequence
!    sfcelev       - true to perform surface obs. elevation check
!    elev_max      - maximum difference between model and obs. height
!    new_sfc_qc    - true to replace NCEP surface QC
!    new_satwnd_qc - true to replace NCEP sat wind QC over ocean
!    rawin_seq     - rawinsonde sequence
!    sfc_seq       - surface sequence
!    acars_seq     - aircraft sequence
!    satwnd_seq    - satellite wind sequence
!    tc_seq        - TC data sequence
!    other_seq     - remaining observation sequence
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_and_parse_input_seq(filename, nx, ny, obs_bdy, siglevel, ptop, &
                                    htop, sfcelev, elev_max, new_sfc_qc, &
                                    new_satwnd_qc, overwrite_time, atime, &
                                    rawin_seq, sfc_seq, acars_seq, satwnd_seq, &
                                    tc_seq, gpsro_seq, other_seq)

use         types_mod, only : r8
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file, &
                                  nc_get_variable
use  time_manager_mod, only : time_type 
use      location_mod, only : location_type, get_location, is_vertical
use  obs_sequence_mod, only : obs_sequence_type, obs_type, init_obs, &
                              get_num_copies, get_num_qc, get_qc_meta_data, &
                              get_first_obs, get_obs_def, copy_obs, get_num_qc, &
                              append_obs_to_seq, get_next_obs, get_qc, set_qc, &
                              destroy_obs_sequence, read_obs_seq, set_obs_def
use       obs_def_mod, only : obs_def_type, get_obs_def_type_of_obs, get_obs_def_location, &
                              set_obs_def_time
use      obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, &
                              RADIOSONDE_SURFACE_ALTIMETER, RADIOSONDE_TEMPERATURE, &
                              RADIOSONDE_SPECIFIC_HUMIDITY, RADIOSONDE_DEWPOINT, &
                              RADIOSONDE_RELATIVE_HUMIDITY, GPSRO_REFRACTIVITY, &
                              AIRCRAFT_U_WIND_COMPONENT, AIRCRAFT_V_WIND_COMPONENT, &
                              AIRCRAFT_TEMPERATURE, AIRCRAFT_SPECIFIC_HUMIDITY, &
                              ACARS_DEWPOINT, ACARS_RELATIVE_HUMIDITY, &
                              ACARS_U_WIND_COMPONENT, ACARS_V_WIND_COMPONENT, &
                              ACARS_TEMPERATURE, ACARS_SPECIFIC_HUMIDITY, &
                              MARINE_SFC_U_WIND_COMPONENT, MARINE_SFC_V_WIND_COMPONENT, &
                              MARINE_SFC_TEMPERATURE, MARINE_SFC_SPECIFIC_HUMIDITY, &
                              MARINE_SFC_RELATIVE_HUMIDITY, MARINE_SFC_DEWPOINT, &
                              LAND_SFC_U_WIND_COMPONENT, LAND_SFC_V_WIND_COMPONENT, &
                              LAND_SFC_TEMPERATURE, LAND_SFC_SPECIFIC_HUMIDITY, &
                              LAND_SFC_RELATIVE_HUMIDITY, LAND_SFC_DEWPOINT, &
                              METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, &
                              METAR_TEMPERATURE_2_METER, METAR_SPECIFIC_HUMIDITY_2_METER, &
                              METAR_DEWPOINT_2_METER, METAR_RELATIVE_HUMIDITY_2_METER, &
                              METAR_ALTIMETER, MARINE_SFC_ALTIMETER, LAND_SFC_ALTIMETER, &
                              SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT, &
                              VORTEX_LAT, VORTEX_LON, VORTEX_PMIN, VORTEX_WMAX
use         model_mod, only : get_domain_info

implicit none

real(r8), parameter                    :: satwnd_qc_ok = 15.0_r8
real(r8), parameter                    :: sfc_qc_ok1   =  9.0_r8
real(r8), parameter                    :: sfc_qc_ok2   = 15.0_r8
real(r8), parameter                    :: new_qc_value =  2.0_r8

character(len=129),      intent(in)    :: filename
real(r8),                intent(in)    :: nx, ny, obs_bdy, ptop, htop, elev_max
logical,                 intent(in)    :: siglevel, sfcelev, new_sfc_qc, &
                                          new_satwnd_qc, overwrite_time
type(time_type),         intent(in)    :: atime
type(obs_sequence_type), intent(inout) :: rawin_seq, sfc_seq, acars_seq, &
                                          satwnd_seq, tc_seq, gpsro_seq, other_seq

character(len=129)    :: qcmeta
integer               :: fid, var_id, okind, dom_id, i, j
logical               :: file_exist, last_obs, input_ncep_qc
real(r8), allocatable :: xland(:,:), qc(:)
real(r8)              :: xyz_loc(3), xloc, yloc

type(location_type)     :: obs_loc
type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: seq
type(obs_type)          :: obs, obs_in, prev_obs

inquire(file = trim(adjustl(filename)), exist = file_exist)
if ( .not. file_exist )  return

call read_obs_seq(filename, 0, 0, 0, seq)

call init_obs(obs,      get_num_copies(seq), get_num_qc(seq))
call init_obs(obs_in,   get_num_copies(seq), get_num_qc(seq))
call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))
allocate(qc(get_num_qc(seq)))

!  read land distribution
allocate(xland(nint(nx),nint(ny)))

fid = nc_open_file_readonly("wrfinput_d01", "read_and_parse_input_seq")
call nc_get_variable(fid, "XLAND", xland)
call nc_close_file(fid, "read_and_parse_input_seq")

input_ncep_qc = .false.
qcmeta = get_qc_meta_data(seq, 1)
if ( trim(adjustl(qcmeta)) == 'NCEP QC index' )  input_ncep_qc = .true.

last_obs = .false.
if ( .not. get_first_obs(seq, obs_in) ) last_obs = .true.

InputObsLoop:  do while ( .not. last_obs ) ! loop over all observations in a sequence

  !  Get the observation information, check if it is in the domain
  call get_obs_def(obs_in, obs_def)
  okind   = get_obs_def_type_of_obs(obs_def)
  obs_loc = get_obs_def_location(obs_def)
  xyz_loc = get_location(obs_loc)
  call get_domain_info(xyz_loc(1),xyz_loc(2),dom_id,xloc,yloc)
  i = nint(xloc);  j = nint(yloc)

  !  check horizontal location
  if ( ((xloc < (obs_bdy+1.0_r8) .or. xloc > (nx-obs_bdy-1.0_r8) .or. &
         yloc < (obs_bdy+1.0_r8) .or. yloc > (ny-obs_bdy-1.0_r8)) .and. &
         (dom_id == 1)) .or. dom_id < 1 ) then

    prev_obs = obs_in
    call get_next_obs(seq, prev_obs, obs_in, last_obs)
    cycle InputObsLoop

  end if

  !  check vertical location
  if ( (is_vertical(obs_loc, "PRESSURE") .and. xyz_loc(3) < ptop) .or. &
       (is_vertical(obs_loc, "HEIGHT")   .and. xyz_loc(3) > htop) ) then

    prev_obs = obs_in
    call get_next_obs(seq, prev_obs, obs_in, last_obs)
    cycle InputObsLoop

  end if

  !  overwrite the observation time with the analysis time if desired
  if ( overwrite_time ) then 
 
    call set_obs_def_time(obs_def, atime)
    call set_obs_def(obs_in, obs_def)
  
  end if

  !  perform platform-specific checks
  select case (okind)

    case ( RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, &
           RADIOSONDE_TEMPERATURE, RADIOSONDE_SPECIFIC_HUMIDITY, &
           RADIOSONDE_DEWPOINT, RADIOSONDE_RELATIVE_HUMIDITY, &
           RADIOSONDE_SURFACE_ALTIMETER)

      if ( rawinsonde_obs_check(obs_loc, okind, siglevel, sfcelev, elev_max) ) then

        call copy_obs(obs, obs_in)
        call append_obs_to_seq(rawin_seq, obs)

      end if

    case ( LAND_SFC_U_WIND_COMPONENT, LAND_SFC_V_WIND_COMPONENT, &
           LAND_SFC_TEMPERATURE, LAND_SFC_SPECIFIC_HUMIDITY, &
           LAND_SFC_RELATIVE_HUMIDITY, LAND_SFC_DEWPOINT,  &
           METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, &
           METAR_TEMPERATURE_2_METER, METAR_SPECIFIC_HUMIDITY_2_METER, &
           METAR_DEWPOINT_2_METER, METAR_RELATIVE_HUMIDITY_2_METER, &
           METAR_ALTIMETER, MARINE_SFC_U_WIND_COMPONENT,  &
           MARINE_SFC_V_WIND_COMPONENT, MARINE_SFC_TEMPERATURE, &
           MARINE_SFC_SPECIFIC_HUMIDITY, MARINE_SFC_DEWPOINT, &
           MARINE_SFC_RELATIVE_HUMIDITY, LAND_SFC_ALTIMETER, MARINE_SFC_ALTIMETER )

      if ( surface_obs_check(sfcelev, elev_max, xyz_loc) ) then

        call copy_obs(obs, obs_in)
        if ( new_sfc_qc .and. okind /= LAND_SFC_ALTIMETER .and. &
             okind /= METAR_ALTIMETER .and. okind /= MARINE_SFC_ALTIMETER ) then

          call get_qc(obs, qc)
          if ( (qc(1) == sfc_qc_ok1 .or. qc(1) == sfc_qc_ok2) .and. input_ncep_qc ) then
            qc(1) = new_qc_value
            call set_qc(obs, qc)
          end if

        end if
        call append_obs_to_seq(sfc_seq, obs)

      endif

    case ( AIRCRAFT_U_WIND_COMPONENT, AIRCRAFT_V_WIND_COMPONENT, &
           AIRCRAFT_TEMPERATURE, AIRCRAFT_SPECIFIC_HUMIDITY, &
           ACARS_RELATIVE_HUMIDITY, ACARS_DEWPOINT, &
           ACARS_U_WIND_COMPONENT, ACARS_V_WIND_COMPONENT, &
           ACARS_TEMPERATURE, ACARS_SPECIFIC_HUMIDITY )

      if ( aircraft_obs_check() ) then

        call copy_obs(obs, obs_in)
        call append_obs_to_seq(acars_seq, obs)

      end if

    case ( SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT )

      if ( sat_wind_obs_check() ) then

        call copy_obs(obs, obs_in)
        if ( new_satwnd_qc ) then

          call get_qc(obs, qc)
          if ( qc(1) == satwnd_qc_ok .and. input_ncep_qc .and. &
                                         xland(i,j) > 1.0_r8 ) then
            qc(1) = new_qc_value
            call set_qc(obs, qc)
          end if

        end if
        call append_obs_to_seq(satwnd_seq, obs)

      endif

    case ( VORTEX_LAT, VORTEX_LON, VORTEX_PMIN, VORTEX_WMAX )

      call copy_obs(obs, obs_in)
      call append_obs_to_seq(tc_seq, obs)

    case ( GPSRO_REFRACTIVITY )

      call copy_obs(obs, obs_in)
      call append_obs_to_seq(gpsro_seq, obs)

    case default

      call copy_obs(obs, obs_in)
      call append_obs_to_seq(other_seq, obs)

  end select

  prev_obs = obs_in
  call get_next_obs(seq, prev_obs, obs_in, last_obs)

end do InputObsLoop
call destroy_obs_sequence(seq)

return
end subroutine read_and_parse_input_seq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   remove_sondes_near_tc - subroutine that removes all rawinsonde
!                           observations within a certain distance of
!                           a TC center.
!
!    obs_seq_tc    - TC observation sequence
!    obs_seq_rawin - rawinsonde observation sequence
!    sonde_radii   - observation removal distance
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine remove_sondes_near_tc(seq_tc, seq_rawin, sonde_radii)

use         types_mod, only : r8, earth_radius
use  obs_sequence_mod, only : obs_sequence_type, init_obs, get_num_copies, &
                              get_num_qc, get_first_obs, get_obs_def, &
                              get_next_obs, delete_obs_from_seq, obs_type, &
                              get_next_obs_from_key, get_obs_key
use       obs_def_mod, only : obs_def_type, get_obs_def_location
use      location_mod, only : location_type, operator(==), get_dist

implicit none

type(obs_sequence_type), intent(in)    :: seq_tc
type(obs_sequence_type), intent(inout) :: seq_rawin
real(r8),                intent(in)    :: sonde_radii

integer :: numtc, n
logical :: last_obs, not_in_list, use_obs, first_obs

type(location_type) :: obs_loc, loctc(20)
type(obs_def_type)  :: obs_def
type(obs_type)      :: obs, prev_obs

write(6,*)  'Removing Sonde Data near TC'
call init_obs(obs,      get_num_copies(seq_rawin), get_num_qc(seq_rawin))
call init_obs(prev_obs, get_num_copies(seq_rawin), get_num_qc(seq_rawin))

last_obs = .false.  ;  numtc = 0
if ( .not. get_first_obs(seq_tc, obs) ) last_obs = .true.

! loop over all TC observations, find locations
do while ( .not. last_obs )

  call get_obs_def(obs, obs_def)
  obs_loc = get_obs_def_location(obs_def)
  not_in_list = .true.
  do n = 1, numtc
    if ( obs_loc == loctc(n) )  not_in_list = .false.
  end do
  if ( not_in_list ) then
    numtc        = numtc + 1
    loctc(numtc) = obs_loc
  end if

  prev_obs = obs
  call get_next_obs(seq_tc, prev_obs, obs, last_obs)

end do

if ( numtc == 0 )  return

last_obs = .false.  ;  first_obs = .true.
if ( .not. get_first_obs(seq_rawin, obs) ) last_obs = .true.
do while ( .not. last_obs )  !  loop over all rawinsonde obs, remove too close to TC

  call get_obs_def(obs, obs_def)
  obs_loc = get_obs_def_location(obs_def)

  use_obs = .true.
  do n = 1, numtc
    if ( (get_dist(obs_loc,loctc(n),2,2,.true.) * earth_radius) <= sonde_radii ) use_obs = .false.
  end do

  if ( use_obs ) then

    prev_obs = obs
    call get_next_obs(seq_rawin, prev_obs, obs, last_obs)
    first_obs = .false.

  else

    if ( first_obs ) then
      call delete_obs_from_seq(seq_rawin, obs)
      if( .not. get_first_obs(seq_rawin, obs) )  return
    else
      call delete_obs_from_seq(seq_rawin, obs)
      call get_next_obs_from_key(seq_rawin, get_obs_key(prev_obs), obs, last_obs)
    end if

  end if

end do

return
end subroutine remove_sondes_near_tc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   sat_wind_obs_check - function that determines whether to include an
!                        satellite wind observation in the sequence.
!                        For now, this function is a placeholder and 
!                        returns true.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function sat_wind_obs_check()

use     types_mod, only : r8

implicit none

logical  :: sat_wind_obs_check

sat_wind_obs_check = .true.

return
end function sat_wind_obs_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   superob_aircraft_data - subroutine that creates superobs of 
!                           aircraft data based on the given
!                           horizontal and vertical intervals.
!
!    seq   - aircraft observation sequence
!    hdist - horizontal interval of superobs
!    vdist - vertical interval of superobs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine superob_aircraft_data(seq, atime, hdist, vdist)

use         types_mod, only : r8, missing_r8, earth_radius
use  time_manager_mod, only : time_type
use      location_mod, only : location_type, get_dist, operator(==), &
                              get_location, VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, init_obs, &
                              get_num_copies, get_num_qc, get_first_obs, &
                              get_next_obs, destroy_obs_sequence, get_qc, &
                              get_num_obs, get_obs_values, get_obs_def, &
                              append_obs_to_seq
use       obs_def_mod, only : obs_def_type, get_obs_def_location, &
                              get_obs_def_type_of_obs, get_obs_def_error_variance, &
                              get_obs_def_time
use      obs_kind_mod, only : AIRCRAFT_U_WIND_COMPONENT, ACARS_U_WIND_COMPONENT, &
                              AIRCRAFT_V_WIND_COMPONENT, ACARS_V_WIND_COMPONENT, &
                              AIRCRAFT_TEMPERATURE, ACARS_TEMPERATURE, &
                              AIRCRAFT_SPECIFIC_HUMIDITY, ACARS_SPECIFIC_HUMIDITY, &
                              ACARS_DEWPOINT, ACARS_RELATIVE_HUMIDITY

implicit none

type(obs_sequence_type), intent(inout) :: seq
type(time_type),         intent(in)    :: atime
real(r8), intent(in)                   :: hdist, vdist

integer             :: num_copies, num_qc, nloc, k, locdex, obs_kind, n, &
                       num_obs, poleward_obs
logical             :: last_obs, close_to_greenwich
real(r8)            :: nuwnd, latu, lonu, preu, uwnd, erru, qcu, nvwnd, latv, &
                       lonv, prev, vwnd, errv, qcv, ntmpk, latt, lont, pret, &
                       tmpk, errt, qct, nqvap, latq, lonq, preq, qvap, errq, &
                       dwpt, errd, qcd, ndwpt, latd, lond, pred, relh, errr, &
                       qcr, nrelh, latr, lonr, prer, qcq, obs_dist, &
                       xyz_loc(3), obs_val(1), qc_val(1), lon_degree_limit
type(location_type) :: obs_loc
type(obs_def_type)  :: obs_def
type(obs_type)      :: obs, prev_obs

type airobs_type

  real(r8)            :: lat, lon, pressure, uwnd, uwnd_err, uwnd_qc, &
                         vwnd, vwnd_err, vwnd_qc, tmpk, tmpk_err, tmpk_qc, &
                         qvap, qvap_err, qvap_qc, dwpt, dwpt_err, dwpt_qc, &
                         relh, relh_err, relh_qc
  type(location_type) :: obs_loc
  type(time_type)     :: time

end type airobs_type

type(airobs_type), allocatable :: airobs(:)

write(6,*) 'Super-Obing Aircraft Data'

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(seq)
num_obs    = get_num_obs(seq)

allocate(airobs(num_obs))
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

last_obs = .false.  ;  nloc = 0   ;   poleward_obs = 0
if ( .not. get_first_obs(seq, obs) )  last_obs = .true.

!  loop over all observations in sequence, add to ACARS observation type
do while ( .not. last_obs )

  call get_obs_values(obs, obs_val, 1)
  call get_qc(obs, qc_val, 1)

  call get_obs_def(obs, obs_def)
  obs_loc  = get_obs_def_location(obs_def)
  obs_kind = get_obs_def_type_of_obs(obs_def)
  xyz_loc  = get_location(obs_loc)

  locdex = -1
  do k = nloc, 1, -1

    if ( obs_loc == airobs(k)%obs_loc ) then
      locdex = k
      exit
    end if

  end do

  if ( locdex < 1 ) then  !  create new observation location type

    ! test if we are within hdist of either pole, and punt for now on those 
    ! obs because we can't accurately average points that wrap the poles.
    ! (count up obs here and print later)
    if (pole_check(xyz_loc(1), xyz_loc(2), hdist)) then
        poleward_obs = poleward_obs + 1
        goto 200
    endif
     
    nloc = nloc + 1
    locdex = nloc

    
    airobs(locdex)%lon = xyz_loc(1)
    airobs(locdex)%lat = xyz_loc(2)
    airobs(locdex)%pressure = xyz_loc(3)
    airobs(locdex)%obs_loc  = obs_loc
    airobs(locdex)%uwnd     = missing_r8
    airobs(locdex)%vwnd     = missing_r8
    airobs(locdex)%tmpk     = missing_r8
    airobs(locdex)%qvap     = missing_r8
    airobs(locdex)%dwpt     = missing_r8
    airobs(locdex)%relh     = missing_r8
    airobs(locdex)%time     = get_obs_def_time(obs_def)

  end if

  !  add observation data to type
  if      ( obs_kind == AIRCRAFT_U_WIND_COMPONENT  .or. obs_kind == ACARS_U_WIND_COMPONENT  ) then

    airobs(locdex)%uwnd     = obs_val(1)
    airobs(locdex)%uwnd_qc  = qc_val(1)
    airobs(locdex)%uwnd_err = get_obs_def_error_variance(obs_def) 

  else if ( obs_kind == AIRCRAFT_V_WIND_COMPONENT  .or. obs_kind == ACARS_V_WIND_COMPONENT  ) then

    airobs(locdex)%vwnd     = obs_val(1)
    airobs(locdex)%vwnd_qc  = qc_val(1)
    airobs(locdex)%vwnd_err = get_obs_def_error_variance(obs_def)

  else if ( obs_kind == AIRCRAFT_TEMPERATURE       .or. obs_kind == ACARS_TEMPERATURE       ) then

    airobs(locdex)%tmpk     = obs_val(1)
    airobs(locdex)%tmpk_qc  = qc_val(1)
    airobs(locdex)%tmpk_err = get_obs_def_error_variance(obs_def)

  else if ( obs_kind == AIRCRAFT_SPECIFIC_HUMIDITY .or. obs_kind == ACARS_SPECIFIC_HUMIDITY )  then

    airobs(locdex)%qvap     = obs_val(1)
    airobs(locdex)%qvap_qc  = qc_val(1)
    airobs(locdex)%qvap_err = get_obs_def_error_variance(obs_def)

  else if ( obs_kind == ACARS_DEWPOINT )  then

    airobs(locdex)%dwpt     = obs_val(1)
    airobs(locdex)%dwpt_qc  = qc_val(1)
    airobs(locdex)%dwpt_err = get_obs_def_error_variance(obs_def)

  else if ( obs_kind == ACARS_RELATIVE_HUMIDITY )  then
  
    airobs(locdex)%relh     = obs_val(1)
    airobs(locdex)%relh_qc  = qc_val(1)
    airobs(locdex)%relh_err = get_obs_def_error_variance(obs_def)


  end if

200 continue   ! come here to skip this obs

  prev_obs = obs
  call get_next_obs(seq, prev_obs, obs, last_obs)

end do

if (poleward_obs > 0) then
   write(6, *) 'WARNING: skipped ', poleward_obs, ' of ', poleward_obs+nloc, ' aircraft obs because'
   write(6, *) 'they were within ', hdist, ' KM of the poles (the superobs distance).'
endif

call destroy_obs_sequence(seq)
call create_new_obs_seq(num_copies, num_qc, num_obs, seq)
call init_obs(obs, num_copies, num_qc)

do k = 1, nloc  !  loop over all observation locations

  nuwnd  = 0.0_r8  ;  latu = 0.0_r8  ;  lonu = 0.0_r8  ;  preu = 0.0_r8
  uwnd   = 0.0_r8  ;  erru = 0.0_r8  ;  qcu  = 0.0_r8
  nvwnd  = 0.0_r8  ;  latv = 0.0_r8  ;  lonv = 0.0_r8  ;  prev = 0.0_r8
  vwnd   = 0.0_r8  ;  errv = 0.0_r8  ;  qcv  = 0.0_r8
  ntmpk  = 0.0_r8  ;  latt = 0.0_r8  ;  lont = 0.0_r8  ;  pret = 0.0_r8
  tmpk   = 0.0_r8  ;  errt = 0.0_r8  ;  qct  = 0.0_r8
  nqvap  = 0.0_r8  ;  latq = 0.0_r8  ;  lonq = 0.0_r8  ;  preq = 0.0_r8
  qvap   = 0.0_r8  ;  errq = 0.0_r8  ;  qcq  = 0.0_r8
  ndwpt  = 0.0_r8  ;  latd = 0.0_r8  ;  lond = 0.0_r8  ;  pred = 0.0_r8
  dwpt   = 0.0_r8  ;  errd = 0.0_r8  ;  qcd  = 0.0_r8
  nrelh  = 0.0_r8  ;  latr = 0.0_r8  ;  lonr = 0.0_r8  ;  prer = 0.0_r8
  relh   = 0.0_r8  ;  errr = 0.0_r8  ;  qcr  = 0.0_r8


  if ( airobs(k)%lat /= missing_r8 ) then  !  create initial superob

    call superob_location_check(airobs(k)%lon, airobs(k)%lat, hdist, &
                                close_to_greenwich, lon_degree_limit)
    if (close_to_greenwich) call wrap_lon(airobs(k)%lon, airobs(k)%lon - lon_degree_limit, &
                                                         airobs(k)%lon + lon_degree_limit)

    if ( airobs(k)%uwnd /= missing_r8 ) then
      nuwnd = nuwnd + 1.0_r8
      latu  = latu  + airobs(k)%lat
      lonu  = lonu  + airobs(k)%lon
      preu  = preu  + airobs(k)%pressure
      uwnd  = uwnd  + airobs(k)%uwnd
      erru  = erru  + airobs(k)%uwnd_err
      qcu   = max(qcu,airobs(k)%uwnd_qc)
    end if

    if ( airobs(k)%vwnd /= missing_r8 ) then
      nvwnd = nvwnd + 1.0_r8
      latv  = latv  + airobs(k)%lat
      lonv  = lonv  + airobs(k)%lon
      prev  = prev  + airobs(k)%pressure
      vwnd  = vwnd  + airobs(k)%vwnd
      errv  = errv  + airobs(k)%vwnd_err
      qcv   = max(qcv,airobs(k)%vwnd_qc)
    end if

    if ( airobs(k)%tmpk /= missing_r8 ) then
      ntmpk = ntmpk + 1.0_r8
      latt  = latt  + airobs(k)%lat
      lont  = lont  + airobs(k)%lon
      pret  = pret  + airobs(k)%pressure
      tmpk  = tmpk  + airobs(k)%tmpk
      errt  = errt  + airobs(k)%tmpk_err
      qct   = max(qct,airobs(k)%tmpk_qc)
    end if

    if ( airobs(k)%qvap /= missing_r8 ) then
      nqvap = nqvap + 1.0_r8
      latq  = latq  + airobs(k)%lat
      lonq  = lonq  + airobs(k)%lon
      preq  = preq  + airobs(k)%pressure
      qvap  = qvap  + airobs(k)%qvap
      errq  = errq  + airobs(k)%qvap_err
      qcq   = max(qcq,airobs(k)%qvap_qc)
    end if

    if ( airobs(k)%dwpt /= missing_r8 ) then
      ndwpt = ndwpt + 1.0_r8
      latd  = latd  + airobs(k)%lat
      lond  = lond  + airobs(k)%lon
      pred  = pred  + airobs(k)%pressure
      dwpt  = dwpt  + airobs(k)%dwpt
      errd  = errd  + airobs(k)%dwpt_err
      qcd   = max(qcd,airobs(k)%dwpt_qc)
    end if

    if ( airobs(k)%relh /= missing_r8 ) then
      nrelh = nrelh + 1.0_r8
      latr  = latr  + airobs(k)%lat
      lonr  = lonr  + airobs(k)%lon
      prer  = prer  + airobs(k)%pressure
      relh  = relh  + airobs(k)%relh
      errr  = errr  + airobs(k)%relh_err
      qcr   = max(qcr,airobs(k)%relh_qc)
    end if

    do n = (k+1), nloc

      if ( airobs(n)%lat /= missing_r8 ) then

        !  add observation to superob if within the horizontal and vertical bounds
        obs_dist = get_dist(airobs(k)%obs_loc, airobs(n)%obs_loc, 2, 2, .true.) * earth_radius
        if ( obs_dist <= hdist .and. abs(airobs(k)%pressure-airobs(n)%pressure) <= vdist ) then

          if (close_to_greenwich) call wrap_lon(airobs(n)%lon, airobs(k)%lon - lon_degree_limit, &
                                                               airobs(k)%lon + lon_degree_limit)

          if ( airobs(n)%uwnd /= missing_r8 ) then
            nuwnd = nuwnd + 1.0_r8
            latu  = latu  + airobs(n)%lat
            lonu  = lonu  + airobs(n)%lon
            preu  = preu  + airobs(n)%pressure
            uwnd  = uwnd  + airobs(n)%uwnd
            erru  = erru  + airobs(n)%uwnd_err
            qcu   = max(qcu,airobs(n)%uwnd_qc)
          end if

          if ( airobs(n)%vwnd /= missing_r8 ) then
            nvwnd = nvwnd + 1.0_r8
            latv  = latv  + airobs(n)%lat
            lonv  = lonv  + airobs(n)%lon
            prev  = prev  + airobs(n)%pressure
            vwnd  = vwnd  + airobs(n)%vwnd
            errv  = errv  + airobs(n)%vwnd_err
            qcv   = max(qcv,airobs(n)%vwnd_qc)
          end if

          if ( airobs(n)%tmpk /= missing_r8 ) then
            ntmpk = ntmpk + 1.0_r8
            latt  = latt  + airobs(n)%lat
            lont  = lont  + airobs(n)%lon
            pret  = pret  + airobs(n)%pressure
            tmpk  = tmpk  + airobs(n)%tmpk
            errt  = errt  + airobs(n)%tmpk_err
            qct   = max(qct,airobs(n)%tmpk_qc)
          end if

          if ( airobs(n)%qvap /= missing_r8 ) then
            nqvap = nqvap + 1.0_r8
            latq  = latq  + airobs(n)%lat
            lonq  = lonq  + airobs(n)%lon
            preq  = preq  + airobs(n)%pressure
            qvap  = qvap  + airobs(n)%qvap
            errq  = errq  + airobs(n)%qvap_err
            qcq   = max(qcq,airobs(n)%qvap_qc)
          end if

          if ( airobs(n)%dwpt /= missing_r8 ) then
            ndwpt = ndwpt + 1.0_r8
            latd  = latd  + airobs(n)%lat
            lond  = lond  + airobs(n)%lon
            pred  = pred  + airobs(n)%pressure
            dwpt  = dwpt  + airobs(n)%dwpt
            errd  = errd  + airobs(n)%dwpt_err
            qcd   = max(qcd,airobs(n)%dwpt_qc)
          end if

          if ( airobs(n)%relh /= missing_r8 ) then
            nrelh = nrelh + 1.0_r8
            latr  = latr  + airobs(n)%lat
            lonr  = lonr  + airobs(n)%lon
            prer  = prer  + airobs(n)%pressure
            relh  = relh  + airobs(n)%relh
            errr  = errr  + airobs(n)%relh_err
            qcr   = max(qcr,airobs(n)%relh_qc)
          end if

          airobs(n)%lat = missing_r8

        end if

      end if

    end do

    if ( nuwnd > 0.0_r8 ) then  !  write zonal wind superob

      latu = latu / nuwnd
      lonu = lonu / nuwnd
      if ( lonu >= 360.0_r8 )  lonu = lonu - 360.0_r8
      preu = preu / nuwnd
      uwnd = uwnd / nuwnd
      erru = erru / nuwnd

      call create_obs_type(latu, lonu, preu, VERTISPRESSURE, uwnd, &
                           ACARS_U_WIND_COMPONENT, erru, qcu, atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

    if ( nvwnd > 0.0_r8 ) then  !  write meridional wind superob

      latv = latv / nvwnd
      lonv = lonv / nvwnd
      if ( lonv >= 360.0_r8 )  lonv = lonv - 360.0_r8
      prev = prev / nvwnd
      vwnd = vwnd / nvwnd
      errv = errv / nvwnd

      call create_obs_type(latv, lonv, prev, VERTISPRESSURE, vwnd, &
                           ACARS_V_WIND_COMPONENT, errv, qcv, atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

    if ( ntmpk > 0.0_r8 ) then  !  write temperature superob

      latt = latt / ntmpk
      lont = lont / ntmpk
      if ( lont >= 360.0_r8 )  lont = lont - 360.0_r8
      pret = pret / ntmpk
      tmpk = tmpk / ntmpk
      errt = errt / ntmpk

      call create_obs_type(latt, lont, pret, VERTISPRESSURE, tmpk, & 
                           ACARS_TEMPERATURE, errt, qct, atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

    if ( nqvap > 0.0_r8 ) then  !  write qvapor superob

      latq = latq / nqvap
      lonq = lonq / nqvap
      if ( lonq >= 360.0_r8 )  lonq = lonq - 360.0_r8
      preq = preq / nqvap
      qvap = qvap / nqvap
      errq = errq / nqvap

      call create_obs_type(latq, lonq, preq, VERTISPRESSURE, qvap, & 
                           ACARS_SPECIFIC_HUMIDITY, errq, qcq, atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

    if ( ndwpt > 0.0_r8 ) then  !  write dewpoint temperature superob

      latd = latd / ndwpt
      lond = lond / ndwpt
      if ( lond >= 360.0_r8 )  lond = lond - 360.0_r8
      pred = pred / ndwpt
      dwpt = dwpt / ndwpt
      errd = errd / ndwpt

      call create_obs_type(latd, lond, pred, VERTISPRESSURE, dwpt, &
                           ACARS_DEWPOINT, errd, qcd, atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

    if ( nrelh > 0.0_r8 ) then  !  write relative humidity superob

      latr = latr / nrelh
      lonr = lonr / nrelh
      if ( lonr >= 360.0_r8 )  lonr = lonr - 360.0_r8
      prer = prer / nrelh
      relh = relh / nrelh
      errr = errr / nrelh

      call create_obs_type(latr, lonr, prer, VERTISPRESSURE, relh, &
                           ACARS_RELATIVE_HUMIDITY, errr, qcr, atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

  end if

end do 

return
end subroutine superob_aircraft_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   superob_sat_wind_data - subroutine that creates superobs of 
!                           satellite wind data based on the given
!                           horizontal and vertical intervals.
!
!    seq   - satellite wind observation sequence
!    hdist - horizontal interval of superobs
!    vdist - vertical interval of superobs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine superob_sat_wind_data(seq, atime, hdist, vdist)

use         types_mod, only : r8, missing_r8, earth_radius
use  time_manager_mod, only : time_type
use      location_mod, only : location_type, get_dist, operator(==), &
                              get_location, VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, init_obs, &
                              get_num_copies, get_num_qc, get_first_obs, &
                              get_next_obs, destroy_obs_sequence, get_qc, &
                              get_num_obs, get_obs_values, get_obs_def, &
                              append_obs_to_seq
use       obs_def_mod, only : obs_def_type, get_obs_def_location, &
                              get_obs_def_type_of_obs, get_obs_def_error_variance, &
                              get_obs_def_time
use      obs_kind_mod, only : SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT

implicit none

type(obs_sequence_type), intent(inout) :: seq
type(time_type),         intent(in)    :: atime
real(r8), intent(in)                   :: hdist, vdist

integer             :: num_copies, num_qc, nloc, k, locdex, obs_kind, n, &
                       num_obs, poleward_obs
logical             :: last_obs, close_to_greenwich
real(r8)            :: nwnd, lat, lon, pres, uwnd, erru, qcu, vwnd, &
                       errv, qcv, obs_dist, xyz_loc(3), obs_val(1), qc_val(1), &
                       lon_degree_limit

type(location_type) :: obs_loc
type(obs_def_type)  :: obs_def
type(obs_type)      :: obs, prev_obs

type satobs_type

  real(r8)            :: lat, lon, pressure, uwnd, uwnd_err, uwnd_qc, &
                         vwnd, vwnd_err, vwnd_qc
  type(location_type) :: obs_loc
  type(time_type)     :: time

end type satobs_type

type(satobs_type), allocatable :: satobs(:)

write(6,*) 'Super-Obing Satellite Wind Data'

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(seq)
num_obs    = get_num_obs(seq)

allocate(satobs(num_obs/2))
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

last_obs = .false.  ;  nloc = 0  ;  poleward_obs = 0
if ( .not. get_first_obs(seq, obs) )  last_obs = .true.

!  loop over satellite winds, create list
do while ( .not. last_obs )

  call get_obs_values(obs, obs_val, 1)
  call get_qc(obs, qc_val, 1)

  call get_obs_def(obs, obs_def)
  obs_loc  = get_obs_def_location(obs_def)
  obs_kind = get_obs_def_type_of_obs(obs_def)
  xyz_loc  = get_location(obs_loc)

  !  determine if observation exists
  locdex = -1
  do k = nloc, 1, -1

    if ( obs_loc == satobs(k)%obs_loc ) then
      locdex = k
      exit
    end if

  end do

  if ( locdex < 1 ) then  !  create new observation type

    ! test if we are within hdist of either pole, and punt for now on those 
    ! obs because we can't accurately average points that wrap the poles.
    ! (hdist is radius, in KM, of region of interest.)
    if (pole_check(xyz_loc(1), xyz_loc(2), hdist)) then
        ! count up obs here and print later
        poleward_obs = poleward_obs + 1
        goto 200
    endif
     
    nloc = nloc + 1
    locdex = nloc

    satobs(locdex)%lon = xyz_loc(1)
    satobs(locdex)%lat = xyz_loc(2)
    satobs(locdex)%pressure = xyz_loc(3)
    satobs(locdex)%obs_loc  = obs_loc
    satobs(locdex)%uwnd     = missing_r8
    satobs(locdex)%vwnd     = missing_r8
    satobs(locdex)%time     = get_obs_def_time(obs_def)

  end if

  !  add observation information
  if ( obs_kind == SAT_U_WIND_COMPONENT ) then

    satobs(locdex)%uwnd     = obs_val(1)
    satobs(locdex)%uwnd_qc  = qc_val(1)
    satobs(locdex)%uwnd_err = get_obs_def_error_variance(obs_def) 

  else if ( obs_kind == SAT_V_WIND_COMPONENT ) then

    satobs(locdex)%vwnd     = obs_val(1)
    satobs(locdex)%vwnd_qc  = qc_val(1)
    satobs(locdex)%vwnd_err = get_obs_def_error_variance(obs_def)

  end if

200 continue   ! come here to skip this obs

  prev_obs = obs
  call get_next_obs(seq, prev_obs, obs, last_obs)

end do

if (poleward_obs > 0) then
   write(6, *) 'WARNING: skipped ', poleward_obs, ' of ', poleward_obs+nloc, ' satwind obs because'
   write(6, *) 'they were within ', hdist, ' KM of the poles (the superobs distance).'
endif

!  create new sequence
call destroy_obs_sequence(seq)
call create_new_obs_seq(num_copies, num_qc, num_obs, seq)
call init_obs(obs, num_copies, num_qc)

do k = 1, nloc  ! loop over all locations

  if ( satobs(k)%uwnd /= missing_r8 .and. satobs(k)%vwnd /= missing_r8 ) then

    call superob_location_check(satobs(k)%lon, satobs(k)%lat, hdist, &
                                close_to_greenwich, lon_degree_limit)
    if (close_to_greenwich) call wrap_lon(satobs(k)%lon, satobs(k)%lon - lon_degree_limit, &
                                                         satobs(k)%lon + lon_degree_limit)

    nwnd = 1.0_r8
    lat  = satobs(k)%lat 
    lon  = satobs(k)%lon
    pres = satobs(k)%pressure
    uwnd = satobs(k)%uwnd
    erru = satobs(k)%uwnd_err
    qcu  = satobs(k)%uwnd_qc
    vwnd = satobs(k)%vwnd
    errv = satobs(k)%vwnd_err
    qcv  = satobs(k)%vwnd_qc

    do n = (k+1), nloc  !  loop over remaining obs

      if ( satobs(n)%uwnd /= missing_r8 .and. satobs(n)%vwnd /= missing_r8 ) then

        !  if observation is within horizontal and vertical interval, add it to obs.
        obs_dist = get_dist(satobs(k)%obs_loc, satobs(n)%obs_loc, 2, 2, .true.) * earth_radius
        if ( obs_dist <= hdist .and. abs(satobs(k)%pressure-satobs(n)%pressure) <= vdist ) then

          if (close_to_greenwich) call wrap_lon(satobs(n)%lon, satobs(k)%lon - lon_degree_limit, &
                                                               satobs(k)%lon + lon_degree_limit)

          nwnd = nwnd + 1.0_r8
          lat  = lat  + satobs(n)%lat
          lon  = lon  + satobs(n)%lon
          pres = pres + satobs(n)%pressure
          uwnd = uwnd + satobs(n)%uwnd
          erru = erru + satobs(n)%uwnd_err
          qcu  = max(qcu,satobs(n)%uwnd_qc)
          vwnd = vwnd + satobs(n)%vwnd
          errv = errv + satobs(n)%vwnd_err
          qcv  = max(qcv,satobs(n)%vwnd_qc)

          satobs(n)%uwnd = missing_r8
          satobs(n)%vwnd = missing_r8

        end if

      end if

    end do

    !  create superobs
    lat  = lat  / nwnd
    lon  = lon  / nwnd
    if ( lon >= 360.0_r8 )  lon = lon - 360.0_r8
    pres = pres / nwnd
    uwnd = uwnd / nwnd
    erru = erru / nwnd
    vwnd = vwnd / nwnd
    errv = errv / nwnd

    !  add to observation sequence
    call create_obs_type(lat, lon, pres, VERTISPRESSURE, uwnd, &
                         SAT_U_WIND_COMPONENT, erru, qcu, atime, obs)
    call append_obs_to_seq(seq, obs)

    call create_obs_type(lat, lon, pres, VERTISPRESSURE, vwnd, &
                         SAT_V_WIND_COMPONENT, errv, qcv, atime, obs)
    call append_obs_to_seq(seq, obs)

  end if

end do 

return
end subroutine superob_sat_wind_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   surface_obs_check - function that determines whether to include an
!                       surface observation in the sequence.
!
!    elev_check - true to check elevation difference
!    elev_max   - maximum difference between model and obs. elevation
!    xyz_loc    - longitude, latitude and elevation array
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function surface_obs_check(elev_check, elev_max, xyz_loc)

use     types_mod, only : r8
use  obs_kind_mod, only : QTY_SURFACE_ELEVATION
use     model_mod, only : model_interpolate
use  location_mod, only : set_location, VERTISSURFACE

implicit none

logical, intent(in)  :: elev_check
real(r8), intent(in) :: xyz_loc(3), elev_max

integer              :: istatus(1)
logical              :: surface_obs_check
real(r8)             :: xmod(1), hsfc(1)

surface_obs_check = .true.

if ( elev_check ) then

  call model_interpolate(dummy_ens, 1, set_location(xyz_loc(1), xyz_loc(2), &
      xyz_loc(3), VERTISSURFACE), QTY_SURFACE_ELEVATION, hsfc, istatus)
  if ( abs(hsfc(1) - xyz_loc(3)) > elev_max ) surface_obs_check = .false.

end if

return
end function surface_obs_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   pole_check - determine if we are within km_dist of either pole.
!                function returns true if so, false if not.
!
!    lon       - longitude in degrees 
!    lat       - latitude in degrees
!    km_dist   - horizontal superob radius in kilometers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function pole_check(lon, lat, km_dist)

use     types_mod, only : r8, earth_radius
use  location_mod, only : location_type, get_dist, set_location, VERTISUNDEF

implicit none

real(r8), intent(in) :: lon, lat, km_dist
logical              :: pole_check

type(location_type) :: thisloc, pole


! create a point at this lon/lat, and at the nearest pole
thisloc = set_location(lon, lat, 0.0_r8, VERTISUNDEF)
if (lat >= 0) then
   pole = set_location(0.0_r8, 90.0_r8, 0.0_r8, VERTISUNDEF)
else
   pole = set_location(0.0_r8, -90.0_r8, 0.0_r8, VERTISUNDEF)
endif

! are we within km_distance of that pole?
if ( get_dist(thisloc, pole, 1, 1, .true.) * earth_radius <= km_dist ) then
   pole_check = .true.
else
   pole_check = .false.
endif

return
end function pole_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   superob_location_check - determine if a point is close to the greenwich
!                            longitude based on the given latitude and distance.
!                            if this location is close enough to longitude=0
!                            we have to take action in order not to average points 
!                            at 0 and 360 and end up with 180 by mistake.  as long
!                            as we treat all points consistently when doing an average
!                            the exact determination of 'close enough' isn't critical.
!                            (the minimum and maximum extents in longitude for a given
!                            point and radius is left as an exercise for the reader.  
!                            hint, they are not along the same latitude circle.)
!
!    lon              - longitude in degrees (input)
!    lat              - latitude in degrees (input)
!    km_dist          - horizontal superob radius in kilometers (input)
!    near_greenwich   - returns true if the given lon/lat is potentially within 
!                       km_dist of longitude 0 (output)
!    lon_degree_limit - number of degrees along a latitude circle that the
!                       km_dist equates to, plus a tolerance (output)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine superob_location_check(lon, lat, km_dist, near_greenwich, lon_degree_limit)

use     types_mod, only : r8, PI, earth_radius, RAD2DEG, DEG2RAD
use  location_mod, only : location_type, get_dist, set_location, VERTISUNDEF

implicit none

real(r8), intent(in)  :: lon, lat, km_dist
logical,  intent(out) :: near_greenwich
real(r8), intent(out) :: lon_degree_limit

real(r8)            :: lat_radius
real(r8), parameter :: fudge_factor = 1.2_r8   ! add a flat 20% 


! the problem with trying to superob in a circle that's specified as a 
! radius of kilometers is that it isn't parallel with longitude lines as 
! you approach the poles.  also when averaging lon values near the prime 
! meridian some values are < 360 and some are > 0 but are close in space.
! simple suggestion for all but the highest latitudes:
!  if dist between lat/lon and lat/0 is < hdist, add 360 to any values >= 0.
!  if the final averaged point >= 360 subtract 360 before setting the location.
! this still isn't good enough as you get closer to the poles; there
! lat/lon averages are a huge pain. hdist could be shorter across the 
! pole and therefore a lon that is 180 degrees away could still be 
! 'close enough' to average in the same superob box.  probably the
! best suggestion for this case is to convert lat/lon into 3d cartesian
! coords, average the x/y/z's separately, and then convert back.
! (no code like this has been implemented here.)
!
! obs_dist_in_km = earth_radius_in_km * dist_in_radians
!  (which is what get_dist(loc1, loc2, 0, 0, .true.) returns)
!
! dist_in_radians = obs_dist_in_km / earth_radius_in_km
!

! figure out how far in degrees km_dist is at this latitude 
! if we traveled along the latitude line.  and since on a sphere
! the actual closest point to a longitude line isn't found by following 
! a latitude circle, add a percentage.  as long as all points are
! treated consistently it doesn't matter if the degree value is exact.

lat_radius = earth_radius * cos(lat*DEG2RAD)
lon_degree_limit = ((km_dist / lat_radius) * RAD2DEG) * fudge_factor

! are we within 'lon_degree_limit' of the greenwich line?
if (lon <= lon_degree_limit .or. (360.0_r8 - lon) <= lon_degree_limit) then
   near_greenwich = .true.
else
   near_greenwich = .false.
endif

return
end subroutine superob_location_check


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   wrap_lon  - update the incoming longitude possibly + 360 degrees if
!               the given limits define a region that crosses long=0.
!               all values should be in units of degrees. 'lon' value
!               should be between westlon and eastlon.
!
!    lon         - longitude to update, returns either unchanged or + 360
!    westlon     - westernmost longitude of region in degrees
!    eastlon     - easternmost longitude of region in degrees
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wrap_lon(lon, westlon, eastlon)

use     types_mod, only : r8

!  uniform way to treat longitude ranges, in degrees, on a globe.
!  adds 360 to the incoming lon if region crosses longitude 0 and
!  given point is east of lon=0.

real(r8), intent(inout) :: lon
real(r8), intent(in)    :: westlon, eastlon

real(r8) :: westl, eastl
real(r8), parameter :: circumf = 360.0_r8

! ensure the region boundaries and target point are between 0 and 360.
! the modulo() function handles negative values ok; mod() does not.
westl = modulo(westlon, circumf)
eastl = modulo(eastlon, circumf)
lon   = modulo(lon,     circumf)

! if the 'region' is the entire globe you can return now.
if (westl == eastl) return

! here's where the magic happens:
! normally the western boundary longitude (westl) has a smaller magnitude than
! the eastern one (eastl).  but westl will be larger than eastl if the region
! of interest crosses the prime meridian. e.g. westl=100, eastl=120 doesn't 
! cross it, westl=340, eastl=10 does.  for regions crossing lon=0, a target lon 
! west of lon=0 should not be changed; a target lon east of lon=0 needs +360 degrees.
! e.g. lon=350 stays unchanged; lon=5 becomes lon=365.

if (westl > eastl .and. lon <= eastl) lon = lon + circumf

return
end subroutine wrap_lon

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

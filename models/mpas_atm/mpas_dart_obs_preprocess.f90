! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program mpas_dart_obs_preprocess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   mpas_dart_obs_preprocess 
! 
!   MPAS-DART utility program that adds observations from supplimental obs sequences.  
!   An actual observation time for each obs can be overwritten as the analysis time,
!   and/or the observations beyond the specified time window can be removed.
!   In addition, this program allows users to do the following functions:
!
!     - remove observations above certain pressure/height levels
!     - remove observations where the model and obs topography are large.
!     - remove significant level rawinsonde data
!     - remove rawinsonde observations near TC core
!     - superob aircraft and satellite wind data
!       - average over observations within each voxel (in 3D).
!       - voxels are defined by mpas grid cells (read from grid_definition_filename 
!         in &model_nml) and the vertical ranges specified in the namelist.
!       - bad observations with qc higher than superob_qc_threshold is not used.
!       - the highest qc and obs error available in each voxel are assigned to
!         the superobed observation.
!       - merge acars and aircraft observations into acars.
!
!     created based on wrf_dart_obs_preprocess by Soyoung Ha, NCAR/MMM (Dec. 2014)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use        types_mod, only : r8, missing_r8, earth_radius, RAD2DEG, DEG2RAD
use    utilities_mod, only : error_handler, E_MSG, find_namelist_in_file, &
                             check_namelist_read, nc_check
use time_manager_mod, only : time_type, operator(>=), operator(<), operator(>), operator(<=), &
                             increment_time, decrement_time, operator(-), operator(+), &
                             set_calendar_type, GREGORIAN, set_time, get_time
use     location_mod, only : location_type, get_location, set_location, get_dist, &
                             VERTISUNDEF, VERTISSURFACE, VERTISPRESSURE, &
                             vert_is_pressure, vert_is_height, operator(==)
use obs_sequence_mod, only : append_obs_to_seq, copy_obs, delete_obs_from_seq, &
                             destroy_obs_sequence, get_first_obs, get_last_obs, &
                             get_next_obs, get_next_obs_from_key, get_num_copies, &
                             get_num_obs, get_num_qc, get_obs_def, get_obs_key, &
                             get_obs_values, get_qc, get_qc_meta_data, init_obs, &
                             insert_obs_in_seq, obs_sequence_type, obs_type, &
                             read_obs_seq, read_obs_seq_header, set_copy_meta_data, &
                             set_obs_def, set_obs_values, set_qc, set_qc_meta_data, &
                             static_init_obs_sequence, write_obs_seq, init_obs_sequence
use      obs_def_mod, only : get_obs_def_error_variance, get_obs_def_location, &
                             get_obs_def_time, get_obs_kind, obs_def_type, &
                             set_obs_def_error_variance, set_obs_def_kind, &
                             set_obs_def_location, set_obs_def_time
use     obs_kind_mod, only : ACARS_DEWPOINT, ACARS_RELATIVE_HUMIDITY, ACARS_SPECIFIC_HUMIDITY, &
                             ACARS_TEMPERATURE, ACARS_U_WIND_COMPONENT, ACARS_V_WIND_COMPONENT, &
                             AIRCRAFT_SPECIFIC_HUMIDITY, AIRCRAFT_TEMPERATURE, AIRCRAFT_U_WIND_COMPONENT, &
                             AIRCRAFT_V_WIND_COMPONENT, GPS_PRECIPITABLE_WATER, GPSRO_REFRACTIVITY, &
                             KIND_SURFACE_ELEVATION, LAND_SFC_ALTIMETER, LAND_SFC_DEWPOINT, &
                             LAND_SFC_RELATIVE_HUMIDITY, LAND_SFC_SPECIFIC_HUMIDITY, LAND_SFC_TEMPERATURE, &
                             LAND_SFC_U_WIND_COMPONENT, LAND_SFC_V_WIND_COMPONENT, MARINE_SFC_ALTIMETER, &
                             MARINE_SFC_DEWPOINT, MARINE_SFC_RELATIVE_HUMIDITY, MARINE_SFC_SPECIFIC_HUMIDITY, &
                             MARINE_SFC_TEMPERATURE, MARINE_SFC_U_WIND_COMPONENT, MARINE_SFC_V_WIND_COMPONENT, &
                             METAR_ALTIMETER, METAR_DEWPOINT_2_METER, METAR_RELATIVE_HUMIDITY_2_METER, &
                             METAR_SPECIFIC_HUMIDITY_2_METER, METAR_TEMPERATURE_2_METER, METAR_U_10_METER_WIND, &
                             METAR_V_10_METER_WIND, PROFILER_U_WIND_COMPONENT, PROFILER_V_WIND_COMPONENT, &
                             RADIOSONDE_DEWPOINT, RADIOSONDE_RELATIVE_HUMIDITY, RADIOSONDE_SPECIFIC_HUMIDITY, &
                             RADIOSONDE_SURFACE_ALTIMETER, RADIOSONDE_TEMPERATURE, RADIOSONDE_U_WIND_COMPONENT, &
                             RADIOSONDE_V_WIND_COMPONENT, SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT, &
                             VORTEX_LAT, VORTEX_LON, VORTEX_PMIN, VORTEX_WMAX
use        model_mod, only : static_init_model, get_grid_dims, get_xland, &
                             model_interpolate, find_closest_cell_center

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
                      gpspw_extra        = 'obs_seq.gpspw',      &
                      trop_cyclone_extra = 'obs_seq.tc'
integer            :: max_num_obs              = 1000000  ! Largest number of obs in one sequence

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
real(r8)           :: aircraft_pres_int        = 2500.0_r8  ! pressure interval for super-ob
integer            :: superob_qc_threshold     = 4          ! reject obs with qc > 4 (applied for both aircraft and satwnd)

!  sat wind specific parameters
logical            :: superob_sat_winds        = .false.    ! super-ob sat wind data
real(r8)           :: sat_wind_pres_int        = 2500.0_r8  ! pressure interval for super-ob
logical            :: overwrite_ncep_satwnd_qc = .false.    ! true to overwrite NCEP QC (see instructions)

!  surface obs. specific parameters
logical            :: overwrite_ncep_sfc_qc    = .false.  ! true to overwrite NCEP QC (see instructions)

! lowest height for GPS REFRACTIVITY (SYHA)
real(r8)           :: gpsro_lowest_meter        = 3000.0      ! remove all obs at lower height

!  overwrite or windowing obs time
logical            :: overwrite_obs_time       = .false.  ! true to overwrite all observation times
logical            :: windowing_obs_time       = .false.  ! true to remove obs beyond the time window
real(r8)           :: windowing_int_hour       = 1.5_r8   ! time window [hr] centered on the analysis time

namelist /mpas_obs_preproc_nml/ file_name_input, file_name_output, max_num_obs,     &
         include_sig_data, superob_aircraft, superob_sat_winds, superob_qc_threshold,   &
         sfc_elevation_check, overwrite_ncep_sfc_qc, overwrite_ncep_satwnd_qc, &
         aircraft_pres_int, sat_wind_pres_int, sfc_elevation_tol,   & 
         obs_pressure_top, obs_height_top, gpsro_lowest_meter, sonde_extra, metar_extra, &
         acars_extra, land_sfc_extra, marine_sfc_extra, sat_wind_extra, profiler_extra, &
         trop_cyclone_extra, gpsro_extra, gpspw_extra, tc_sonde_radii, overwrite_obs_time, &
         windowing_obs_time, windowing_int_hour

!----------------------------------------------------------------------
! Declare other variables
!----------------------------------------------------------------------

character(len=129)      :: obs_seq_read_format
character(len=80)       :: name

integer                 :: io, iunit, fid, var_id, obs_seq_file_id, num_copies, &
                           num_qc, num_obs, max_obs_seq, gday, gsec

logical                 :: file_exist, pre_I_format

type(obs_sequence_type) :: seq_all, seq_rawin, seq_sfc, seq_acars, seq_satwnd, &
                           seq_prof, seq_tc, seq_gpsro, seq_other, seq_gpspw

type(time_type)         :: anal_time

integer :: nCells        = -1  ! Total number of cells making up the grid
integer :: nVertices     = -1  ! Unique points in grid that are corners of cells
integer :: nEdges        = -1  ! Straight lines between vertices making up cells
integer :: nVertLevels   = -1  ! Vertical levels; count of vert cell centers
integer :: vertexDegree  = -1  ! Max number of cells/edges that touch any vertex
integer :: nSoilLevels   = -1  ! Number of soil layers

integer :: dimid, ncid, VarID
real(r8), allocatable :: xland(:)   ! indicator for land (1.0_r8) or ocean (> 1.0_r8)

print*,'Enter target assimilation time (gregorian day, second): '
read*, gday,gsec
call set_calendar_type(GREGORIAN)
anal_time = set_time(gsec, gday)

call find_namelist_in_file("input.nml", "mpas_obs_preproc_nml", iunit)
read(iunit, nml = mpas_obs_preproc_nml, iostat = io)
call check_namelist_read(iunit, io, "mpas_obs_preproc_nml")

call static_init_obs_sequence()
call static_init_model()

call get_grid_dims(nCells, nVertices, nEdges, nVertLevels, vertexDegree, nSoilLevels)

allocate(xland(nCells))
call get_xland(nCells,xland)
!print*,'xland: ',minval(xland),maxval(xland)

!  if obs_seq file exists, read in the data, otherwise, create a blank one.
inquire(file = trim(adjustl(file_name_input)), exist = file_exist)

if ( file_exist ) then

  print*,'file_name_input: ',trim(file_name_input)
  call read_obs_seq_header(file_name_input, num_copies, num_qc, num_obs, max_obs_seq, &
        obs_seq_file_id, obs_seq_read_format, pre_I_format, close_the_file = .true.)
 
  if( max_obs_seq < max_num_obs ) max_obs_seq = max_num_obs

else

  num_copies = 1  ;  num_qc = 1  ;  max_obs_seq = max_num_obs
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
call create_new_obs_seq(num_copies, num_qc, max_obs_seq, seq_gpspw)

!  read input obs_seq file, divide into platforms
call read_and_parse_input_seq(file_name_input, xland,                    &
include_sig_data, obs_pressure_top, obs_height_top, sfc_elevation_check, &
sfc_elevation_tol, overwrite_ncep_sfc_qc, overwrite_ncep_satwnd_qc,      &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour,   &
seq_rawin, seq_sfc, seq_acars, seq_satwnd, seq_tc, seq_gpsro,            &
seq_gpspw, seq_other)

!  add supplimental rawinsonde observations from file
call add_supplimental_obs(sonde_extra, seq_rawin, max_obs_seq, &
RADIOSONDE_U_WIND_COMPONENT, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  add supplimental ACARS observations from file
call add_supplimental_obs(acars_extra, seq_acars, max_obs_seq, &
ACARS_U_WIND_COMPONENT, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  add supplimental marine observations from file
call add_supplimental_obs(marine_sfc_extra, seq_sfc, max_obs_seq, &
MARINE_SFC_U_WIND_COMPONENT, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  add supplimental land surface observations from file
call add_supplimental_obs(land_sfc_extra, seq_sfc, max_obs_seq, &
LAND_SFC_U_WIND_COMPONENT, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  add supplimental metar observations from file
call add_supplimental_obs(metar_extra, seq_sfc, max_obs_seq, &
METAR_U_10_METER_WIND, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  add supplimental satellite wind observations from file
call add_supplimental_obs(sat_wind_extra, seq_satwnd, max_obs_seq, &
SAT_U_WIND_COMPONENT, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  add supplimental profiler observations from file
call add_supplimental_obs(profiler_extra, seq_prof, max_obs_seq, &
PROFILER_U_WIND_COMPONENT, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  add supplimental GPSRO observations from file
call add_supplimental_obs(gpsro_extra, seq_gpsro, max_obs_seq, &
GPSRO_REFRACTIVITY, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  add supplimental GPSPW observations from file
call add_supplimental_obs(gpspw_extra, seq_gpspw, max_obs_seq, &
GPS_PRECIPITABLE_WATER, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  add supplimental tropical cyclone vortex observations from file
call add_supplimental_obs(trop_cyclone_extra, seq_tc, max_obs_seq, &
VORTEX_LAT, include_sig_data, &
obs_pressure_top, obs_height_top, gpsro_lowest_meter, sfc_elevation_check, sfc_elevation_tol, &
overwrite_obs_time, anal_time, windowing_obs_time, windowing_int_hour)

!  remove all sonde observations within radius of TC if desired
if ( tc_sonde_radii > 0.0_r8 ) call remove_sondes_near_tc(seq_tc, & 
                                               seq_rawin, tc_sonde_radii)

!  super-ob ACARS data
if ( superob_aircraft ) call superob_aircraft_data(seq_acars, nCells, anal_time, &
                             aircraft_pres_int, superob_qc_threshold, obs_pressure_top)

!  super-ob satellite wind data
if ( superob_sat_winds ) call superob_sat_wind_data(seq_satwnd, nCells, anal_time, &
                              sat_wind_pres_int, superob_qc_threshold, obs_pressure_top)

print*, 'Number of obs read:'
print*, 'num_rawin:  ', get_num_obs(seq_rawin)
print*, 'num_sfc:    ', get_num_obs(seq_sfc)
print*, 'num_acars:  ', get_num_obs(seq_acars)
print*, 'num_satwnd: ', get_num_obs(seq_satwnd)
print*, 'num_prof:   ', get_num_obs(seq_prof)
print*, 'num_gpsro:  ', get_num_obs(seq_gpsro)
print*, 'num_gpspw:  ', get_num_obs(seq_gpspw)
print*, 'num_tc:     ', get_num_obs(seq_tc)
print*, 'num_other:  ', get_num_obs(seq_other)

max_obs_seq = get_num_obs(seq_tc)     + get_num_obs(seq_rawin) + &
              get_num_obs(seq_sfc)    + get_num_obs(seq_acars) + &
              get_num_obs(seq_satwnd) + get_num_obs(seq_prof)  + &
              get_num_obs(seq_gpsro)  + get_num_obs(seq_gpspw) + &
              get_num_obs(seq_other)
print*, 'num_total:  ', max_obs_seq

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

call build_master_sequence(seq_gpspw, seq_all)
call destroy_obs_sequence(seq_gpspw)

call build_master_sequence(seq_satwnd, seq_all)
call destroy_obs_sequence(seq_satwnd)

call build_master_sequence(seq_prof, seq_all)
call destroy_obs_sequence(seq_prof)

call build_master_sequence(seq_other, seq_all)
call destroy_obs_sequence(seq_other)

write(6,*) 'Total number of observations after superobing:', get_num_obs(seq_all)
write(6,*) ''

!  write the observation sequence to file
call write_obs_seq(seq_all, file_name_output)
call destroy_obs_sequence(seq_all)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   aircraft_obs_check - function that determines whether to include an
!                     aircraft observation in the sequence.  For now,
!                     this function is a placeholder and returns true.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function aircraft_obs_check()

logical  :: aircraft_obs_check

aircraft_obs_check = .true.

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
!    siglevel    - true to include sonde significant level data
!    ptop        - lowest pressure to include in sequence
!    htop        - highest height level to include in sequence
!    sfcelev     - true to perform surface obs. elevation check
!    elev_max    - maximum difference between model and obs. height
!    overwrite_time - if true, replace actual observation time with atime
!    atime       - analysis time, for windowing and overwriting obs times
!    obs_window  - if true, exclude obs earlier or later than window interval
!    window_hours - hours for time window, obs more than +/- away discarded
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_supplimental_obs(filename, obs_seq, max_obs_seq, plat_kind, &
                                 siglevel, ptop, htop, hbot, sfcelev, elev_max,  &
                                 overwrite_time, atime, obs_window, window_hours)

character(len=129),      intent(in)    :: filename
type(obs_sequence_type), intent(inout) :: obs_seq
integer,                 intent(in)    :: max_obs_seq, plat_kind
logical,                 intent(in)    :: siglevel, sfcelev, overwrite_time
real(r8),                intent(in)    :: ptop, htop, hbot, elev_max
type(time_type),         intent(in)    :: atime
logical,                 intent(in)    :: obs_window
real(r8),                intent(in)    :: window_hours

integer  :: nloc, okind, dom_id
integer  :: gsec, gday, dsec, bday, bsec, eday, esec, num_excluded_bytime
logical  :: file_exist, last_obs, pass_checks, first_obs
real(r8) :: llv_loc(3)

type(location_type)     :: obs_loc_list(max_obs_seq), obs_loc
type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: supp_obs_seq
type(obs_type)          :: obs_in, prev_obsi, prev_obso, obs
type(time_type)         :: obs_time, prev_time
type(time_type)         :: window_min, window_max

inquire(file = trim(adjustl(filename)), exist = file_exist)
if ( .not. file_exist )  return

write(6,*) ''
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
  case (GPS_PRECIPITABLE_WATER)
    write(6,*) 'Adding Supplimental GPS PW Data'

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

! windowing obs - don't compute these things if not going to use them
if ( obs_window ) then
  dsec = nint(window_hours * 3600.)
  window_min = decrement_time(atime, dsec)
  window_max = increment_time(atime, dsec)
  num_excluded_bytime    = 0   ! total number of obs beyond the time window
end if

ObsLoop:  do while ( .not. last_obs ) ! loop over all observations in a sequence

  !  read data from observation
  call get_obs_def(obs_in, obs_def)
  okind   = get_obs_kind(obs_def)
  obs_loc = get_obs_def_location(obs_def)
  llv_loc = get_location(obs_loc)
  obs_time = get_obs_def_time(obs_def)

  !  check if the observation is within vertical bounds of domain
  if ( (vert_is_pressure(obs_loc) .and. llv_loc(3) < ptop) .or. &
       (vert_is_height(obs_loc)   .and. llv_loc(3) > htop) ) then

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

  if ( obs_window ) then
    if ( obs_time <= window_min .or. obs_time > window_max ) then

      prev_obsi = obs_in
      call get_next_obs(supp_obs_seq, prev_obsi, obs_in, last_obs)
      num_excluded_bytime = num_excluded_bytime + 1
      cycle ObsLoop

    end if
  end if

  ! perform platform-specific checks
  select case (plat_kind)

    case (RADIOSONDE_U_WIND_COMPONENT)
      pass_checks = rawinsonde_obs_check(obs_loc, okind, siglevel, &
                                                sfcelev, elev_max)
    case (ACARS_U_WIND_COMPONENT)
      pass_checks = aircraft_obs_check()
    case (MARINE_SFC_U_WIND_COMPONENT)
      pass_checks = surface_obs_check(sfcelev, elev_max, llv_loc)
    case (LAND_SFC_U_WIND_COMPONENT)
      pass_checks = surface_obs_check(sfcelev, elev_max, llv_loc)
    case (METAR_U_10_METER_WIND)
      pass_checks = surface_obs_check(sfcelev, elev_max, llv_loc)
    case (SAT_U_WIND_COMPONENT)
      pass_checks = sat_wind_obs_check()
    case (GPSRO_REFRACTIVITY)
      pass_checks = minimum_height_check(hbot, llv_loc)

    case default
      pass_checks = .true.

  end select

  if ( pass_checks ) then

    call copy_obs(obs, obs_in)
    call get_obs_def(obs, obs_def)
    obs_time = get_obs_def_time(obs_def)

    !  overwrite the observation time with the analysis time if desired
    if ( overwrite_time ) then
         call set_obs_def_time(obs_def, atime)
         call set_obs_def(obs, obs_def)
    end if

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

if ( obs_window ) &
print*, 'Number of obs outside the time window in this supplimental_obs:',num_excluded_bytime

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

integer, intent(in)           :: okind, vcord
real(r8), intent(in)          :: lat, lon, vloc, obsv, oerr, qc
type(time_type), intent(in)   :: otime
type(obs_type), intent(inout) :: obs

real(r8)              :: obs_val(1), qc_val(1)
type(obs_def_type)    :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vloc, vcord))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, otime)
call set_obs_def_error_variance(obs_def, oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

end subroutine create_obs_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   isManLevel - function that returns a logical true if the input 
!                pressure level is a mandatory rawinsonde level.
!
!    plevel - pressure level to check (Pa)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function isManLevel(plevel)
real(r8), intent(in) :: plevel

integer, parameter :: nman = 16
integer :: kk
logical :: isManLevel
real (r8) raw_man_levels(nman) &
     / 100000.0_r8, 92500.0_r8, 85000.0_r8, 70000.0_r8, 50000.0_r8, 40000.0_r8, &
        30000.0_r8, 25000.0_r8, 20000.0_r8, 15000.0_r8, 10000.0_r8,  7000.0_r8, &
         5000.0_r8,  3000.0_r8,  2000.0_r8,  1000.0_r8 /

isManLevel = .false.
do kk = 1, nman
  if ( plevel == raw_man_levels(kk) ) then
    isManLevel = .true.
    return 
  end if
end do

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

integer, intent(in)             :: nloc
type(location_type), intent(in) :: obsloc, obsloc_list(nloc)
logical :: original_observation

real(r8), parameter :: dist_epsilon = 0.00001_r8
integer :: n

original_observation = .true.

do n = 1, nloc

  if ( get_dist(obsloc, obsloc_list(n), 1, 1, .true.) <= dist_epsilon ) then
    original_observation = .false.
    return
  end if

end do

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

type(location_type), intent(in) :: obs_loc
integer, intent(in)             :: obs_kind
logical, intent(in)             :: siglevel, elev_check
real(r8), intent(in)            :: elev_max
logical  :: rawinsonde_obs_check

integer  :: istatus
real(r8) :: llv_loc(3), xmod(1), hsfc

rawinsonde_obs_check = .true.
llv_loc = get_location(obs_loc)

if ( obs_kind /= RADIOSONDE_SURFACE_ALTIMETER ) then

  !  check if vertical level is mandatory level
  if ( (.not. siglevel) .and. (.not. isManLevel(llv_loc(3))) ) then
    rawinsonde_obs_check = .false.
    return
  end if

else

  !  perform elevation check for altimeter
  if ( elev_check ) then

    call model_interpolate(xmod, obs_loc, KIND_SURFACE_ELEVATION, hsfc, istatus)
    if ( abs(hsfc - llv_loc(3)) > elev_max ) rawinsonde_obs_check = .false.

  end if

end if

end function rawinsonde_obs_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   read_and_parse_input_seq - subroutine that reads a generic
!                              observation sequence and divides the
!                              obs into sequences for each platform.
!
!    filename      - name of input obs sequence
!    landmask      - land/ocean mask, dimensioned (nCells), 1=land,2=water
!    siglevel      - true to include sonde significant level data
!    ptop          - lowest pressure to include in sequence
!    htop          - highest height level to include in sequence
!    sfcelev       - true to perform surface obs. elevation check
!    elev_max      - maximum difference between model and obs. height
!    new_sfc_qc    - true to replace NCEP surface QC
!    new_satwnd_qc - true to replace NCEP sat wind QC over ocean
!    overwrite_time - if true, replace actual observation time with atime
!    atime       - analysis time, for windowing and overwriting obs times
!    obs_window  - if true, exclude obs earlier or later than window interval
!    window_hours - hours for time window, obs more than +/- away discarded
!    rawin_seq     - rawinsonde sequence
!    sfc_seq       - surface sequence
!    acars_seq     - aircraft sequence
!    satwnd_seq    - satellite wind sequence
!    tc_seq        - TC data sequence
!    gpspw_seq     - total precipitable water from GPS observations
!    other_seq     - remaining observation sequence
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_and_parse_input_seq(filename, landmask, siglevel, ptop,        &
                                    htop, sfcelev, elev_max, new_sfc_qc,       &
                                    new_satwnd_qc, overwrite_time, atime,      &
                                    obs_window, window_hours,                  &
                                    rawin_seq, sfc_seq, acars_seq, satwnd_seq, &
                                    tc_seq, gpsro_seq, gpspw_seq, other_seq)

character(len=129),      intent(in)    :: filename
real(r8),                intent(in)    :: landmask(:)
real(r8),                intent(in)    :: ptop, htop, elev_max
logical,                 intent(in)    :: siglevel, sfcelev, new_sfc_qc, &
                                          new_satwnd_qc, overwrite_time
logical,                 intent(in)    :: obs_window
real(r8),                intent(in)    :: window_hours
type(time_type),         intent(in)    :: atime
type(obs_sequence_type), intent(inout) :: rawin_seq, sfc_seq, acars_seq, gpspw_seq, &
                                          satwnd_seq, tc_seq, gpsro_seq, other_seq

real(r8), parameter :: satwnd_qc_ok = 15.0_r8
real(r8), parameter :: sfc_qc_ok1   =  9.0_r8
real(r8), parameter :: sfc_qc_ok2   = 15.0_r8
real(r8), parameter :: new_qc_value =  2.0_r8

character(len=129)    :: qcmeta
integer               :: fid, var_id, okind, dom_id, cellid, dsec
integer               :: bsec, bday, esec, eday, num_excluded_bytime
logical               :: file_exist, last_obs, input_ncep_qc
real(r8), allocatable :: qc(:)
real(r8)              :: llv_loc(3)

type(location_type)     :: obs_loc
type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: seq
type(obs_type)          :: obs, obs_in, prev_obs
type(time_type)         :: window_min, window_max, obs_time

inquire(file = trim(adjustl(filename)), exist = file_exist)
if ( .not. file_exist )  return

call read_obs_seq(filename, 0, 0, 0, seq)
write(6,*) ''
write(6,*) 'Total of ',get_num_obs(seq),' observations in ',trim(filename)

call init_obs(obs,      get_num_copies(seq), get_num_qc(seq))
call init_obs(obs_in,   get_num_copies(seq), get_num_qc(seq))
call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))
allocate(qc(get_num_qc(seq)))

input_ncep_qc = .false.
qcmeta = get_qc_meta_data(seq, 1)
if ( trim(adjustl(qcmeta)) == 'NCEP QC index' )  input_ncep_qc = .true.

last_obs = .false.
if ( .not. get_first_obs(seq, obs_in) ) last_obs = .true.


if ( obs_window ) then
     write(6,*) 'Time windowing obs within +- ',window_hours,' hours.'
     dsec = nint(window_hours * 3600.)
     window_min = decrement_time(atime, dsec)
     window_max = increment_time(atime, dsec)
     call get_time(window_min,bsec,bday)
     call get_time(window_max,esec,eday)
     print*, 'including obs after ',bday,bsec,' up to and including',eday,esec
     num_excluded_bytime = 0   ! total number of obs beyond the time window
end if


InputObsLoop:  do while ( .not. last_obs ) ! loop over all observations in a sequence

  !  Get the observation information, check if it is in the domain
  call get_obs_def(obs_in, obs_def)
  okind   = get_obs_kind(obs_def)
  obs_loc = get_obs_def_location(obs_def)
  llv_loc = get_location(obs_loc)
  cellid  = find_closest_cell_center(llv_loc(2), llv_loc(1))
  obs_time = get_obs_def_time(obs_def)

  !  check vertical location
  if ( (vert_is_pressure(obs_loc) .and. llv_loc(3) < ptop) .or. &
       (vert_is_height(obs_loc)   .and. llv_loc(3) > htop) ) then

    prev_obs = obs_in
    call get_next_obs(seq, prev_obs, obs_in, last_obs)
    cycle InputObsLoop

  end if

  ! to prevent obs from being included in more than a single assimilation window
  ! discard obs that are on the earliest boundary but keep obs that are on the
  ! latest boundary.  this matches what dart does with assimilation windows; they
  ! run from start+1 second to end.
  if ( obs_window ) then
    if (obs_time <= window_min .or. obs_time > window_max ) then
         prev_obs = obs_in
         call get_next_obs(seq, prev_obs, obs_in, last_obs)
         num_excluded_bytime = num_excluded_bytime + 1
         cycle InputObsLoop
     end if
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

      if ( surface_obs_check(sfcelev, elev_max, llv_loc) ) then

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
               landmask(cellid) > 1.0_r8 ) then    ! LAND MASK (1 FOR LAND, 2 FOR WATER) - not over land
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

    case ( GPS_PRECIPITABLE_WATER )

      call copy_obs(obs, obs_in)
      call append_obs_to_seq(gpspw_seq, obs)

    case default

      call copy_obs(obs, obs_in)
      call append_obs_to_seq(other_seq, obs)

  end select

  prev_obs = obs_in
  call get_next_obs(seq, prev_obs, obs_in, last_obs)

end do InputObsLoop
call destroy_obs_sequence(seq)
if ( obs_window ) &
print*, 'Number of obs outside the time window in the input file:',num_excluded_bytime

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
logical  :: sat_wind_obs_check

sat_wind_obs_check = .true.

end function sat_wind_obs_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   superob_aircraft_data - subroutine that creates superobs of 
!                           aircraft data based on the given
!                           horizontal and vertical intervals.
!
!    seq   - aircraft observation sequence
!    vdist - vertical interval of superobs
!    ptop  - lowest pressure to include in sequence
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine superob_aircraft_data(seq, ncell, atime, vdist, iqc_thres, ptop)

type(obs_sequence_type), intent(inout) :: seq
type(time_type),         intent(in)    :: atime
integer,  intent(in)                   :: ncell
integer,  intent(in)                   :: iqc_thres   ! Superob obs w/ qc < iqc_thres only.
real(r8), intent(in)                   :: vdist, ptop

character(len=256)  :: string
integer             :: icell
integer             :: num_copies, num_qc, nloc, k, locdex, obs_kind, n, &
                       num_obs, poleward_obs
logical             :: last_obs
real(r8)            :: llv_loc(3), obs_val(1), qc_val(1)
real(r8),allocatable:: nuwnd(:,:),latu(:,:),lonu(:,:),preu(:,:),uwnd(:,:),erru(:,:),qcu(:,:),&
                       nvwnd(:,:),latv(:,:),lonv(:,:),prev(:,:),vwnd(:,:),errv(:,:),qcv(:,:),& 
                       ntmpk(:,:),latt(:,:),lont(:,:),pret(:,:),tmpk(:,:),errt(:,:),qct(:,:),&
                       nqvap(:,:),latq(:,:),lonq(:,:),preq(:,:),qvap(:,:),errq(:,:),qcq(:,:),&
                       ndwpt(:,:),latd(:,:),lond(:,:),pred(:,:),dwpt(:,:),errd(:,:),qcd(:,:),&
                       nrelh(:,:),latr(:,:),lonr(:,:),prer(:,:),relh(:,:),errr(:,:),qcr(:,:)

integer             :: nlev, ik
real(r8)            :: ps, pt, dp
real(r8),allocatable:: plevs(:)

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

!-----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
!-----------------------------------------------------------------------

write(6,*)
write(6,*) 'Super-Obing Aircraft Data over ',ncell,' cells.'

! Vertical layers for superobing up to obs_pressure_top.
! plevs is defined at the midpoint between two adjent levels.
! We do not use plevs, but define it here just in case one wants to 
! check the levels to which aircraft obs were superobed.
ps = 100000.0_r8    ! Pa
pt = ptop
dp = vdist * 2
nlev = nint((ps - pt)/dp) + 1

allocate(plevs(nlev))
do k = 1, nlev
   plevs(k) = ps - dp * ( k - 1 )
enddo
write(6,*) 'Super-Obing to',nlev, ' pressure levels (Pa) vertically.'
write(6,'(5f10.1)') plevs(1:nlev-1)

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(seq)
num_obs    = get_num_obs(seq)

write(6,*) 'Super-Obing',num_obs,' Aircraft Data'

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
  obs_kind = get_obs_kind(obs_def)
  llv_loc  = get_location(obs_loc)

  locdex = -1
  do k = nloc, 1, -1

    if ( obs_loc == airobs(k)%obs_loc ) then
      locdex = k
      exit
    end if

  end do

  if ( locdex < 1 ) then  !  create new observation location type

    ! test if we are within the cell of either pole, and punt for now on those 
    ! obs because we can't accurately average points that wrap the poles.
    ! (count up obs here and print later)
    if (pole_check(llv_loc(1), llv_loc(2))) then
        poleward_obs = poleward_obs + 1
        goto 200
    endif
     
    nloc = nloc + 1
    locdex = nloc
    
    airobs(locdex)%lon      = llv_loc(1)
    airobs(locdex)%lat      = llv_loc(2)
    airobs(locdex)%pressure = llv_loc(3)
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

   if(qc_val(1).lt.iqc_thres) then
    airobs(locdex)%uwnd     = obs_val(1)
    airobs(locdex)%uwnd_qc  = qc_val(1)
    airobs(locdex)%uwnd_err = get_obs_def_error_variance(obs_def) 
   endif

  else if ( obs_kind == AIRCRAFT_V_WIND_COMPONENT  .or. obs_kind == ACARS_V_WIND_COMPONENT  ) then

   if(qc_val(1).lt.iqc_thres) then
    airobs(locdex)%vwnd     = obs_val(1)
    airobs(locdex)%vwnd_qc  = qc_val(1)
    airobs(locdex)%vwnd_err = get_obs_def_error_variance(obs_def)
   endif

  else if ( obs_kind == AIRCRAFT_TEMPERATURE       .or. obs_kind == ACARS_TEMPERATURE       ) then

   if(qc_val(1).lt.iqc_thres) then
    airobs(locdex)%tmpk     = obs_val(1)
    airobs(locdex)%tmpk_qc  = qc_val(1)
    airobs(locdex)%tmpk_err = get_obs_def_error_variance(obs_def)
   endif

  else if ( obs_kind == AIRCRAFT_SPECIFIC_HUMIDITY .or. obs_kind == ACARS_SPECIFIC_HUMIDITY )  then

   if(qc_val(1).lt.iqc_thres) then
    airobs(locdex)%qvap     = obs_val(1)
    airobs(locdex)%qvap_qc  = qc_val(1)
    airobs(locdex)%qvap_err = get_obs_def_error_variance(obs_def)
   endif

  else if ( obs_kind == ACARS_DEWPOINT )  then

   if(qc_val(1).lt.iqc_thres) then
    airobs(locdex)%dwpt     = obs_val(1)
    airobs(locdex)%dwpt_qc  = qc_val(1)
    airobs(locdex)%dwpt_err = get_obs_def_error_variance(obs_def)
   endif

  else if ( obs_kind == ACARS_RELATIVE_HUMIDITY )  then
  
   if(qc_val(1).lt.iqc_thres) then
    airobs(locdex)%relh     = obs_val(1)
    airobs(locdex)%relh_qc  = qc_val(1)
    airobs(locdex)%relh_err = get_obs_def_error_variance(obs_def)
   endif

  end if

200 continue   ! come here to skip this obs

  prev_obs = obs
  call get_next_obs(seq, prev_obs, obs, last_obs)

end do         !do while ( .not. last_obs )

if (poleward_obs > 0) then
   write(6, *) 'WARNING: skipped ', poleward_obs, ' of ', poleward_obs+nloc, ' aircraft obs because'
   write(6, *) 'they were within the grid cell at the poles.'
   !write(6, *) 'they were within ', hdist, ' KM of the poles (the superobs distance).'
endif

call destroy_obs_sequence(seq)
call create_new_obs_seq(num_copies, num_qc, num_obs, seq)
call init_obs(obs, num_copies, num_qc)

! Allocation and initialization
allocate(nuwnd(ncell,nlev)); allocate(nvwnd(ncell,nlev)); allocate(ntmpk(ncell,nlev));
allocate(nqvap(ncell,nlev)); allocate(ndwpt(ncell,nlev)); allocate(nrelh(ncell,nlev));
allocate(latu(ncell,nlev));  allocate(latv(ncell,nlev));  allocate(latt(ncell,nlev));
allocate(latq(ncell,nlev));  allocate(latd(ncell,nlev));  allocate(latr(ncell,nlev));
allocate(lonu(ncell,nlev));  allocate(lonv(ncell,nlev));  allocate(lont(ncell,nlev));
allocate(lonq(ncell,nlev));  allocate(lond(ncell,nlev));  allocate(lonr(ncell,nlev));
allocate(preu(ncell,nlev));  allocate(prev(ncell,nlev));  allocate(pret(ncell,nlev));
allocate(preq(ncell,nlev));  allocate(pred(ncell,nlev));  allocate(prer(ncell,nlev));
allocate(uwnd(ncell,nlev));  allocate(vwnd(ncell,nlev));  allocate(tmpk(ncell,nlev));
allocate(qvap(ncell,nlev));  allocate(dwpt(ncell,nlev));  allocate(relh(ncell,nlev));
allocate(erru(ncell,nlev));  allocate(errv(ncell,nlev));  allocate(errt(ncell,nlev));
allocate(errq(ncell,nlev));  allocate(errd(ncell,nlev));  allocate(errr(ncell,nlev));
allocate( qcu(ncell,nlev));  allocate( qcv(ncell,nlev));  allocate( qct(ncell,nlev));
allocate( qcq(ncell,nlev));  allocate( qcd(ncell,nlev));  allocate( qcr(ncell,nlev));

nuwnd=0.0_r8; latu=0.0_r8; lonu=0.0_r8; preu=0.0_r8; uwnd=0.0_r8; erru=0.0_r8; qcu=0.0_r8
nvwnd=0.0_r8; latv=0.0_r8; lonv=0.0_r8; prev=0.0_r8; vwnd=0.0_r8; errv=0.0_r8; qcv=0.0_r8
ntmpk=0.0_r8; latt=0.0_r8; lont=0.0_r8; pret=0.0_r8; tmpk=0.0_r8; errt=0.0_r8; qct=0.0_r8
nqvap=0.0_r8; latq=0.0_r8; lonq=0.0_r8; preq=0.0_r8; qvap=0.0_r8; errq=0.0_r8; qcq=0.0_r8
ndwpt=0.0_r8; latd=0.0_r8; lond=0.0_r8; pred=0.0_r8; dwpt=0.0_r8; errd=0.0_r8; qcd=0.0_r8
nrelh=0.0_r8; latr=0.0_r8; lonr=0.0_r8; prer=0.0_r8; relh=0.0_r8; errr=0.0_r8; qcr=0.0_r8


! Assign obs into each bin [ncell, nlev] for superobing
do k = 1, nloc  !  loop over all observation locations

   icell=0;    ik=0
   icell = find_closest_cell_center(airobs(k)%lat, airobs(k)%lon)
   if(icell < 1) then
      write(string,*) 'Cannot find any cell for this obs at ', airobs(k)%lat, airobs(k)%lon
      call error_handler(E_MSG,'superob_aircraft_data', source, revision, revdate) 
   endif
   if(airobs(k)%pressure > ps) then
      ik = 1
   else   
      ik = nint((ps - airobs(k)%pressure)/dp) + 1
   endif

   if ( airobs(k)%uwnd /= missing_r8 ) then
      nuwnd(icell,ik) = nuwnd(icell,ik) + 1.0_r8
       latu(icell,ik) = latu(icell,ik)  + airobs(k)%lat
       lonu(icell,ik) = lonu(icell,ik)  + airobs(k)%lon
       preu(icell,ik) = preu(icell,ik)  + airobs(k)%pressure
       uwnd(icell,ik) = uwnd(icell,ik)  + airobs(k)%uwnd
       erru(icell,ik) = max(erru(icell,ik),airobs(k)%uwnd_err)
      !erru(icell,ik) = erru(icell,ik)  + airobs(k)%uwnd_err
        qcu(icell,ik) = max(qcu(icell,ik),airobs(k)%uwnd_qc)
   end if

   if ( airobs(k)%vwnd /= missing_r8 ) then
      nvwnd(icell,ik) = nvwnd(icell,ik) + 1.0_r8
       latv(icell,ik) = latv(icell,ik)  + airobs(k)%lat
       lonv(icell,ik) = lonv(icell,ik)  + airobs(k)%lon
       prev(icell,ik) = prev(icell,ik)  + airobs(k)%pressure
       vwnd(icell,ik) = vwnd(icell,ik)  + airobs(k)%vwnd
       errv(icell,ik) = max(errv(icell,ik),airobs(k)%vwnd_err)
      !errv(icell,ik) = errv(icell,ik)  + airobs(k)%vwnd_err
        qcv(icell,ik) = max(qcv(icell,ik),airobs(k)%vwnd_qc)
   end if

   if ( airobs(k)%tmpk /= missing_r8 ) then
      ntmpk(icell,ik) = ntmpk(icell,ik) + 1.0_r8
       latt(icell,ik) = latt(icell,ik)  + airobs(k)%lat
       lont(icell,ik) = lont(icell,ik)  + airobs(k)%lon
       pret(icell,ik) = pret(icell,ik)  + airobs(k)%pressure
       tmpk(icell,ik) = tmpk(icell,ik)  + airobs(k)%tmpk
       errt(icell,ik) = max(errt(icell,ik),airobs(k)%tmpk_err)
      !errt(icell,ik) = errt(icell,ik)  + airobs(k)%tmpk_err
        qct(icell,ik) = max(qct(icell,ik),airobs(k)%tmpk_qc)
   end if

   if ( airobs(k)%qvap /= missing_r8 ) then
      nqvap(icell,ik) = nqvap(icell,ik) + 1.0_r8
       latq(icell,ik) = latq(icell,ik)  + airobs(k)%lat
       lonq(icell,ik) = lonq(icell,ik)  + airobs(k)%lon
       preq(icell,ik) = preq(icell,ik)  + airobs(k)%pressure
       qvap(icell,ik) = qvap(icell,ik)  + airobs(k)%qvap
       errq(icell,ik) = max(errq(icell,ik),airobs(k)%qvap_err)
      !errq(icell,ik) = errq(icell,ik)  + airobs(k)%qvap_err
       qcq(icell,ik)  = max(qcq(icell,ik),airobs(k)%qvap_qc)
   end if

   if ( airobs(k)%dwpt /= missing_r8 ) then
      ndwpt(icell,ik) = ndwpt(icell,ik) + 1.0_r8
       latd(icell,ik) = latd(icell,ik)  + airobs(k)%lat
       lond(icell,ik) = lond(icell,ik)  + airobs(k)%lon
       pred(icell,ik) = pred(icell,ik)  + airobs(k)%pressure
       dwpt(icell,ik) = dwpt(icell,ik)  + airobs(k)%dwpt
       errd(icell,ik) = max(errd(icell,ik),airobs(k)%dwpt_err)
      !errd(icell,ik) = errd(icell,ik)  + airobs(k)%dwpt_err
        qcd(icell,ik) = max(qcd(icell,ik),airobs(k)%dwpt_qc)
   end if

   if ( airobs(k)%relh /= missing_r8 ) then
      nrelh(icell,ik) = nrelh(icell,ik) + 1.0_r8
       latr(icell,ik) = latr(icell,ik)  + airobs(k)%lat
       lonr(icell,ik) = lonr(icell,ik)  + airobs(k)%lon
       prer(icell,ik) = prer(icell,ik)  + airobs(k)%pressure
       relh(icell,ik) = relh(icell,ik)  + airobs(k)%relh
       errr(icell,ik) = max(errr(icell,ik),airobs(k)%relh_err)
      !errr(icell,ik) = errr(icell,ik)  + airobs(k)%relh_err
        qcr(icell,ik) = max(qcr(icell,ik),airobs(k)%relh_qc)
   end if

enddo !k = 1, nloc  !  loop over all observation locations
deallocate(airobs)

do k = 1, nlev   ! loop over all vertical levels

do n = 1, ncell  ! loop over all grid cells

   if ( nuwnd(n,k) > 0.0_r8 ) then !  write zonal wind superob

        latu(n,k) = latu(n,k) / nuwnd(n,k) 
        lonu(n,k) = lonu(n,k) / nuwnd(n,k)
        if ( lonu(n,k) >= 360.0_r8 )  lonu(n,k) = lonu(n,k) - 360.0_r8
        preu(n,k) = preu(n,k) / nuwnd(n,k)
        uwnd(n,k) = uwnd(n,k) / nuwnd(n,k)
       !erru(n,k) = erru(n,k) / nuwnd(n,k)

        call create_obs_type(latu(n,k), lonu(n,k), preu(n,k), VERTISPRESSURE, uwnd(n,k), &
                             ACARS_U_WIND_COMPONENT, erru(n,k), qcu(n,k), atime, obs)
        call append_obs_to_seq(seq, obs)

   end if

   if ( nvwnd(n,k) > 0.0_r8 ) then  !  write meridional wind superob

      latv(n,k) = latv(n,k) / nvwnd(n,k)
      lonv(n,k) = lonv(n,k) / nvwnd(n,k)
      if ( lonv(n,k) >= 360.0_r8 )  lonv(n,k) = lonv(n,k) - 360.0_r8
      prev(n,k) = prev(n,k) / nvwnd(n,k)
      vwnd(n,k) = vwnd(n,k) / nvwnd(n,k)
     !errv(n,k) = errv(n,k) / nvwnd(n,k)

      call create_obs_type(latv(n,k), lonv(n,k), prev(n,k), VERTISPRESSURE, vwnd(n,k), &
                           ACARS_V_WIND_COMPONENT, errv(n,k), qcv(n,k), atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

    if ( ntmpk(n,k) > 0.0_r8 ) then  !  write temperature superob

      latt(n,k) = latt(n,k) / ntmpk(n,k)
      lont(n,k) = lont(n,k) / ntmpk(n,k)
      if ( lont(n,k) >= 360.0_r8 )  lont(n,k) = lont(n,k) - 360.0_r8
      pret(n,k) = pret(n,k) / ntmpk(n,k)
      tmpk(n,k) = tmpk(n,k) / ntmpk(n,k)
     !errt(n,k) = errt(n,k) / ntmpk(n,k)

      call create_obs_type(latt(n,k), lont(n,k), pret(n,k), VERTISPRESSURE, tmpk(n,k), & 
                           ACARS_TEMPERATURE, errt(n,k), qct(n,k), atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

    if ( nqvap(n,k) > 0.0_r8 ) then  !  write qvapor superob

      latq(n,k) = latq(n,k) / nqvap(n,k)
      lonq(n,k) = lonq(n,k) / nqvap(n,k)
      if ( lonq(n,k) >= 360.0_r8 )  lonq(n,k) = lonq(n,k) - 360.0_r8
      preq(n,k) = preq(n,k) / nqvap(n,k)
      qvap(n,k) = qvap(n,k) / nqvap(n,k)
     !errq(n,k) = errq(n,k) / nqvap(n,k)

      call create_obs_type(latq(n,k), lonq(n,k), preq(n,k), VERTISPRESSURE, qvap(n,k), & 
                           ACARS_SPECIFIC_HUMIDITY, errq(n,k), qcq(n,k), atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

    if ( ndwpt(n,k) > 0.0_r8 ) then  !  write dewpoint temperature superob

      latd(n,k) = latd(n,k) / ndwpt(n,k)
      lond(n,k) = lond(n,k) / ndwpt(n,k)
      if ( lond(n,k) >= 360.0_r8 )  lond(n,k) = lond(n,k) - 360.0_r8
      pred(n,k) = pred(n,k) / ndwpt(n,k)
      dwpt(n,k) = dwpt(n,k) / ndwpt(n,k)
     !errd(n,k) = errd(n,k) / ndwpt(n,k)

      call create_obs_type(latd(n,k), lond(n,k), pred(n,k), VERTISPRESSURE, dwpt(n,k), &
                           ACARS_DEWPOINT, errd(n,k), qcd(n,k), atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

    if ( nrelh(n,k) > 0.0_r8 ) then  !  write relative humidity superob

      latr(n,k) = latr(n,k) / nrelh(n,k)
      lonr(n,k) = lonr(n,k) / nrelh(n,k)
      if ( lonr(n,k) >= 360.0_r8 )  lonr(n,k) = lonr(n,k) - 360.0_r8
      prer(n,k) = prer(n,k) / nrelh(n,k)
      relh(n,k) = relh(n,k) / nrelh(n,k)
     !errr(n,k) = errr(n,k) / nrelh(n,k)

      call create_obs_type(latr(n,k), lonr(n,k), prer(n,k), VERTISPRESSURE, relh(n,k), &
                           ACARS_RELATIVE_HUMIDITY, errr(n,k), qcr(n,k), atime, obs)
      call append_obs_to_seq(seq, obs)

    end if

enddo 

enddo !k = 1, nlev   ! loop over all vertical levels

deallocate(plevs)
deallocate(nuwnd); deallocate(nvwnd); deallocate(ntmpk);
deallocate(nqvap); deallocate(ndwpt); deallocate(nrelh);
deallocate(latu);  deallocate(latv);  deallocate(latt);
deallocate(latq);  deallocate(latd);  deallocate(latr);
deallocate(lonu);  deallocate(lonv);  deallocate(lont);
deallocate(lonq);  deallocate(lond);  deallocate(lonr);
deallocate(preu);  deallocate(prev);  deallocate(pret);
deallocate(preq);  deallocate(pred);  deallocate(prer);
deallocate(uwnd);  deallocate(vwnd);  deallocate(tmpk);
deallocate(qvap);  deallocate(dwpt);  deallocate(relh);
deallocate(erru);  deallocate(errv);  deallocate(errt);
deallocate(errq);  deallocate(errd);  deallocate(errr);
deallocate( qcu);  deallocate( qcv);  deallocate( qct);
deallocate( qcq);  deallocate( qcd);  deallocate( qcr);

end subroutine superob_aircraft_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   superob_sat_wind_data - subroutine that creates superobs of 
!                           satellite wind data based on the given
!                           horizontal and vertical intervals.
!
!    seq   - satellite wind observation sequence
!    vdist - vertical interval of superobs
!    ptop  - lowest pressure to include in sequence
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine superob_sat_wind_data(seq, ncell, atime, vdist, iqc_thres, ptop)


type(obs_sequence_type), intent(inout) :: seq
type(time_type),         intent(in)    :: atime
integer,  intent(in)                   :: ncell, iqc_thres
real(r8), intent(in)                   :: vdist, ptop

character(len=256)  :: string
integer             :: icell
integer             :: num_copies, num_qc, nloc, k, locdex, obs_kind, n, &
                       num_obs, poleward_obs
logical             :: last_obs
real(r8)            :: llv_loc(3), obs_val(1), qc_val(1)

real(r8),allocatable :: nwnd(:,:), lat(:,:), lon(:,:), pres(:,:), &
                        uwnd(:,:), erru(:,:), qcu(:,:), vwnd(:,:), errv(:,:), qcv(:,:)
real(r8),allocatable :: plevs(:)
real(r8)             :: ps, pt, dp
integer              :: nlev, ik

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

!-----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
!-----------------------------------------------------------------------

write(6,*)
write(6,*) 'Super-Obing Satellite Wind Data over ',ncell,' cells.'

! Vertical layers for superobing up to obs_pressure_top.
! plevs is defined at the midpoint between two adjent levels.
! We do not use plevs, but define it here just in case one wants to 
! check the levels to which aircraft obs were superobed.
ps = 100000.0_r8    ! Pa
pt = ptop
dp = vdist * 2
nlev = nint((ps - pt)/dp) + 1

allocate(plevs(nlev))
do k = 1, nlev
   plevs(k) = ps - dp * ( k - 1 )
enddo

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(seq)
num_obs    = get_num_obs(seq)

write(6,*) 'Super-Obing',num_obs,' Satellite Wind Data'
write(6,*)

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
  obs_kind = get_obs_kind(obs_def)
  llv_loc  = get_location(obs_loc)

  !  determine if observation exists
  locdex = -1
  do k = nloc, 1, -1

    if ( obs_loc == satobs(k)%obs_loc ) then
      locdex = k
      exit
    end if

  end do

  if ( locdex < 1 ) then  !  create new observation type

    ! test if we are within the cell of either pole, and punt for now on those 
    ! obs because we can't accurately average points that wrap the poles.
    ! (hdist is radius, in KM, of region of interest.)
    if (pole_check(llv_loc(1), llv_loc(2))) then
        ! count up obs here and print later
        poleward_obs = poleward_obs + 1
        goto 200
    endif
     
    nloc = nloc + 1
    locdex = nloc

    satobs(locdex)%lon      = llv_loc(1)
    satobs(locdex)%lat      = llv_loc(2)
    satobs(locdex)%pressure = llv_loc(3)
    satobs(locdex)%obs_loc  = obs_loc
    satobs(locdex)%uwnd     = missing_r8
    satobs(locdex)%vwnd     = missing_r8
    satobs(locdex)%time     = get_obs_def_time(obs_def)

  end if

  !  add observation information
  if ( obs_kind == SAT_U_WIND_COMPONENT ) then

   if(qc_val(1).lt.iqc_thres) then
    satobs(locdex)%uwnd     = obs_val(1)
    satobs(locdex)%uwnd_qc  = qc_val(1)
    satobs(locdex)%uwnd_err = get_obs_def_error_variance(obs_def) 
   endif

  else if ( obs_kind == SAT_V_WIND_COMPONENT ) then

   if(qc_val(1).lt.iqc_thres) then
    satobs(locdex)%vwnd     = obs_val(1)
    satobs(locdex)%vwnd_qc  = qc_val(1)
    satobs(locdex)%vwnd_err = get_obs_def_error_variance(obs_def)
   endif

  end if

200 continue   ! come here to skip this obs

  prev_obs = obs
  call get_next_obs(seq, prev_obs, obs, last_obs)

end do


if (poleward_obs > 0) then
   write(6, *) 'WARNING: skipped ', poleward_obs, ' of ', poleward_obs+nloc, ' satwind obs because'
   write(6, *) 'they were within the cell of the poles (the superobs distance).'
   !write(6, *) 'they were within ', hdist, ' cells of the poles (the superobs distance).'
endif

!  create new sequence
call destroy_obs_sequence(seq)
call create_new_obs_seq(num_copies, num_qc, num_obs, seq)
call init_obs(obs, num_copies, num_qc)

! Allocation and initialization
allocate(nwnd(ncell,nlev)); allocate(pres(ncell,nlev))
allocate( lat(ncell,nlev)); allocate( lon(ncell,nlev))
allocate(uwnd(ncell,nlev)); allocate(vwnd(ncell,nlev))
allocate(erru(ncell,nlev)); allocate(errv(ncell,nlev))
allocate( qcu(ncell,nlev)); allocate( qcv(ncell,nlev))

nwnd=0.0_r8; pres=0.0_r8;  lat=0.0_r8;  lon=0.0_r8
uwnd=0.0_r8; vwnd=0.0_r8; erru=0.0_r8; errv=0.0_r8
 qcu=0.0_r8;  qcv=0.0_r8

icell=0;    ik=0

! Assign obs into each bin [ncell, nlev] for superobing
do k = 1, nloc  ! loop over all locations

   if ( satobs(k)%uwnd /= missing_r8 .and. satobs(k)%vwnd /= missing_r8 ) then

        icell = find_closest_cell_center(satobs(k)%lat, satobs(k)%lon)
        if(icell < 1) then
           write(string,*) 'Cannot find any cell for this obs at ',&
                            satobs(k)%lat,satobs(k)%lon  
           call error_handler(E_MSG,'superob_sat_wind_data', source, revision, revdate) 
        endif
        if(satobs(k)%pressure > ps) then
           ik = 1
        else
           ik = nint((ps - satobs(k)%pressure)/dp) + 1
        endif

        nwnd(icell,ik) = nwnd(icell,ik) + 1.0_r8
         lat(icell,ik) =  lat(icell,ik) + satobs(k)%lat 
         lon(icell,ik) =  lon(icell,ik) + satobs(k)%lon 
        pres(icell,ik) = pres(icell,ik) + satobs(k)%pressure 
        uwnd(icell,ik) = uwnd(icell,ik) + satobs(k)%uwnd 
        vwnd(icell,ik) = vwnd(icell,ik) + satobs(k)%vwnd 
        erru(icell,ik) = erru(icell,ik) + satobs(k)%uwnd_err
        errv(icell,ik) = errv(icell,ik) + satobs(k)%vwnd_err
         qcu(icell,ik) = max(qcu(icell,ik),satobs(k)%uwnd_qc)
         qcv(icell,ik) = max(qcv(icell,ik),satobs(k)%vwnd_qc)

   end if

end do    !  do k = 1, nloc  ! loop over all locations

! Superob in each bin [ncell, nlev]
do k = 1, nlev  ! loop over all locations

do n = 1, ncell ! loop over all grid cells

   if( nwnd(n,k) > 0.0_r8 ) then      ! superob

     ! create superobs
       lat(n,k)  = lat(n,k)  / nwnd(n,k)
       lon(n,k)  = lon(n,k)  / nwnd(n,k)
       if ( lon(n,k) >= 360.0_r8 )  lon(n,k) = lon(n,k) - 360.0_r8
       pres(n,k) = pres(n,k) / nwnd(n,k)
       uwnd(n,k) = uwnd(n,k) / nwnd(n,k)
       erru(n,k) = erru(n,k) / nwnd(n,k)
       vwnd(n,k) = vwnd(n,k) / nwnd(n,k)
       errv(n,k) = errv(n,k) / nwnd(n,k)

       ! NCEP satwnd over land
       if(qcu(n,k).eq.9 .or. qcu(n,k).eq.15) then
          qcu(n,k) = 0
       endif
       if(qcv(n,k).eq.9 .or. qcv(n,k).eq.15) then
          qcv(n,k) = 0
       endif

     ! add to observation sequence
       call create_obs_type(lat(n,k), lon(n,k), pres(n,k), VERTISPRESSURE, uwnd(n,k), &
                            SAT_U_WIND_COMPONENT, erru(n,k), qcu(n,k), atime, obs)
       call append_obs_to_seq(seq, obs)

       call create_obs_type(lat(n,k), lon(n,k), pres(n,k), VERTISPRESSURE, vwnd(n,k), &
                            SAT_V_WIND_COMPONENT, errv(n,k), qcv(n,k), atime, obs)
       call append_obs_to_seq(seq, obs)

   endif     !( nwnd(n,k) > 0.0_r8 ) then      ! superob

end do       ! do n = 1, ncell

end do       ! do k = 1, nlev

deallocate(plevs)
deallocate(nwnd); deallocate(pres) 
deallocate(uwnd); deallocate(vwnd)
deallocate(erru); deallocate(errv)
deallocate(lat); deallocate(lon)
deallocate(qcu); deallocate(qcv)

end subroutine superob_sat_wind_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   minimum_height_check - function that determines whether to include an
!                          observation based on whether the height is
!                          above the given minimum
!
!    min_height - lowest accepted height (in meters)
!    llv_loc    - longitude, latitude and elevation array
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function minimum_height_check(min_height, llv_loc)

real(r8), intent(in) :: llv_loc(3), min_height
logical              :: minimum_height_check

minimum_height_check = .true.

if (llv_loc(3) < min_height) then

   minimum_height_check = .false.

endif

end function minimum_height_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   surface_obs_check - function that determines whether to include an
!                       surface observation in the sequence.
!
!    elev_check - true to check elevation difference
!    elev_max   - maximum difference between model and obs. elevation
!    llv_loc    - longitude, latitude and elevation array
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function surface_obs_check(elev_check, elev_max, llv_loc)

logical, intent(in)  :: elev_check
real(r8), intent(in) :: llv_loc(3), elev_max

integer              :: istatus
logical              :: surface_obs_check
real(r8)             :: xmod(1), hsfc

surface_obs_check = .true.

if ( elev_check ) then

  call model_interpolate(xmod, set_location(llv_loc(1), llv_loc(2), &
      llv_loc(3), VERTISSURFACE), KIND_SURFACE_ELEVATION, hsfc, istatus)
  if ( abs(hsfc - llv_loc(3)) > elev_max ) surface_obs_check = .false.

end if

end function surface_obs_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   pole_check - determine if we are within the cell enclosing either pole.
!                function returns true if so, false if not.
!
!    lon       - longitude in degrees 
!    lat       - latitude in degrees
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function pole_check(lon, lat)

real(r8), intent(in) :: lon, lat
logical              :: pole_check
logical, save        :: first = .true.
integer, save        :: north_pole, south_pole

integer              :: cellid, poleid

if(first) then
   first = .false.
   north_pole = find_closest_cell_center( 90.0_r8, 0.0_r8)
   south_pole = find_closest_cell_center(-90.0_r8, 0.0_r8)
endif

! create a point at this lon/lat, and at the nearest pole
cellid  = find_closest_cell_center(lat, lon)

! are we within the cell at that pole?
! FIXME: For now, we check if the obs is located within the cell at the pole.
! Later on, we may use hdist to account for multiple cells like in superobing.
! Can we make hist an optional argument?
if ( cellid .eq. north_pole .or. cellid .eq. south_pole ) then
   pole_check = .true.
else
   pole_check = .false.
endif

end function pole_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   superob_location_check - determine if a point is located within the predefined
!                            number of cells.
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

real(r8), intent(in)  :: lon, lat, km_dist
logical,  intent(out) :: near_greenwich
real(r8), intent(out) :: lon_degree_limit

real(r8)            :: lat_radius
real(r8), parameter :: fudge_factor = 1.2_r8   ! add a flat 20% 



lat_radius = earth_radius * cos(lat*DEG2RAD)
lon_degree_limit = ((km_dist / lat_radius) * RAD2DEG) * fudge_factor

! are we within 'lon_degree_limit' of the greenwich line?
if (lon <= lon_degree_limit .or. (360.0_r8 - lon) <= lon_degree_limit) then
   near_greenwich = .true.
else
   near_greenwich = .false.
endif

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

end subroutine wrap_lon

end program


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
PROGRAM gts_to_dart

use         types_mod, only : r8, missing_r8, missing_data, DEG2RAD
use     utilities_mod, only : close_file, initialize_utilities, &
                              register_module, logfileunit, E_MSG, finalize_utilities, &
                              error_handler, find_namelist_in_file, check_namelist_read
use  obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, insert_obs_in_seq, &
                              set_copy_meta_data, set_qc_meta_data, write_obs_seq, assignment(=), &
                              init_obs, static_init_obs_sequence, set_obs_def, set_obs_values, set_qc
use       obs_def_mod, only : set_obs_def_location, set_obs_def_error_variance, &
                              set_obs_def_type_of_obs, set_obs_def_time, set_obs_def_key, &
                              obs_def_type
use      obs_kind_mod, only : SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT, &
                              QKSWND_U_WIND_COMPONENT, QKSWND_V_WIND_COMPONENT, &
                              RADIOSONDE_U_WIND_COMPONENT, RADIOSONDE_V_WIND_COMPONENT, &
                              RADIOSONDE_TEMPERATURE, RADIOSONDE_SPECIFIC_HUMIDITY, &
                              RADIOSONDE_DEWPOINT, METAR_DEWPOINT_2_METER, &
                              METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, METAR_TEMPERATURE_2_METER, &
                              METAR_SPECIFIC_HUMIDITY_2_METER, METAR_SURFACE_PRESSURE, METAR_POT_TEMP_2_METER, &
                              BUOY_U_WIND_COMPONENT, BUOY_V_WIND_COMPONENT, BUOY_SURFACE_PRESSURE, &
                              BUOY_TEMPERATURE, BUOY_DEWPOINT, &
                              SHIP_U_WIND_COMPONENT, SHIP_V_WIND_COMPONENT, SHIP_SURFACE_PRESSURE, &
                              SHIP_TEMPERATURE, SHIP_DEWPOINT, &
                              SYNOP_U_WIND_COMPONENT, SYNOP_V_WIND_COMPONENT, SYNOP_SURFACE_PRESSURE, &
                              SYNOP_TEMPERATURE, SYNOP_DEWPOINT, &
                              METAR_U_10_METER_WIND, METAR_V_10_METER_WIND, METAR_SURFACE_PRESSURE, &
                              METAR_TEMPERATURE_2_METER, &
                              AIREP_U_WIND_COMPONENT, AIREP_V_WIND_COMPONENT, AIREP_PRESSURE, &
                              AIREP_TEMPERATURE, AIREP_DEWPOINT, &
                              AMDAR_U_WIND_COMPONENT, AMDAR_V_WIND_COMPONENT, AMDAR_PRESSURE, &
                              AMDAR_TEMPERATURE, AMDAR_DEWPOINT, &
                              PILOT_U_WIND_COMPONENT, PILOT_V_WIND_COMPONENT, PILOT_PRESSURE, &
                              PILOT_TEMPERATURE, PILOT_DEWPOINT, &
                              PROFILER_U_WIND_COMPONENT, PROFILER_V_WIND_COMPONENT, PROFILER_PRESSURE, &
                              BOGUS_U_WIND_COMPONENT, BOGUS_V_WIND_COMPONENT, BOGUS_PRESSURE, &
                              BOGUS_TEMPERATURE, BOGUS_DEWPOINT, &
                              GPSRO_REFRACTIVITY
use      location_mod, only : location_type, set_location, VERTISSURFACE, VERTISPRESSURE, VERTISHEIGHT
use  time_manager_mod, only : time_type, set_date, set_calendar_type, GREGORIAN

use DA_Constants
use DA_Define_Structures
use module_obs

use gts_dart_mod

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: obs
type(obs_def_type)      :: obs_def
type(location_type)     :: location
type(time_type)         :: time

integer           :: iunit, iost, io
integer           :: year, month, day, hours, minutes, seconds
integer           :: i, ii, which_vert, gpskey
integer           :: num_obs, num_obs_last, num_copies, num_qc, max_num_obs=8000000

!real(r8)          :: lat,lon,pressure, height
!real(r8)          :: windspd, winddir, uu, vv, tk, td, rh, qv

character(len=80) :: dummy

logical           :: valid

!-----------------------------------------------------------------------------
! Namelist with default values
!
character(len = 129) :: gts_file        = 'gts_obsout.dat', &
                        obs_seq_out_file_name = 'obs_seq.out'
integer              :: calendar_type         = GREGORIAN

integer              :: gts_qc_threshold = -1     
                        ! < 0 : keep all obs (see qc_flag in da_control in WRFVAR)
                        ! > 0 : keep obs with 0 <= qc_flag <= gts_qc_threshold


logical  ::  Use_SynopObs   = .TRUE., &
             Use_ShipsObs   = .TRUE., &
             Use_MetarObs   = .TRUE., &
             Use_BuoysObs   = .TRUE., &
             Use_PilotObs   = .TRUE., &
             Use_SoundObs   = .TRUE., &
             Use_SatemObs   = .TRUE., &
             Use_SatobObs   = .TRUE., &
             Use_AirepObs   = .TRUE., &
             Use_AmdarObs   = .TRUE., &
             Use_GpspwObs   = .TRUE., &
             Use_SsmiRetrievalObs = .TRUE., &
             Use_SsmiTbObs  = .TRUE., &
             Use_Ssmt1Obs   = .TRUE., &
             Use_Ssmt2Obs   = .TRUE., &
             Use_ProflObs   = .TRUE., &
             Use_QscatObs   = .TRUE., &
             Use_BogusObs   = .TRUE., &
             Use_gpsrefobs   = .TRUE.
logical  ::  dropsonde_only  = .FALSE.

integer  ::  num_thin_satob = 50, num_thin_qscat = 100

namelist /gts_to_dart_nml/ gts_file, obs_seq_out_file_name, &
                           gts_qc_threshold, &
                           Use_SynopObs, &
                           Use_ShipsObs, &
                           Use_MetarObs, &
                           Use_BuoysObs, &
                           Use_PilotObs, &
                           Use_SoundObs, &
                           Use_SatemObs, &
                           Use_SatobObs, &
                           Use_AirepObs, &
                           Use_AmdarObs, &
                           Use_GpspwObs, &
                           Use_SsmiRetrievalObs, &
                           Use_SsmiTbObs, &
                           Use_Ssmt1Obs, &
                           Use_Ssmt2Obs, &
                           Use_QscatObs, &
                           Use_ProflObs, &
                           Use_gpsrefobs, &
                           Use_Bogusobs, &
                           dropsonde_only, &
                           num_thin_satob, &
                           num_thin_qscat 

!------------------------------------------------------------------------------

!print*,'Input WRFVAR 2.1 GTS OBS filename: '
!read(*,'(a)') gts_file
!print*,'Output DART OBS Sequence filename: '
!read(*,'(a)') obs_seq_out_file_name

call initialize_utilities('gts_to_dart')
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "gts_to_dart_nml", iunit)
read(iunit, nml = gts_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "gts_to_dart_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'gts_to_dart','gts_to_dart_nml values are',' ',' ',' ')
write(logfileunit, nml=gts_to_dart_nml)
write(     *     , nml=gts_to_dart_nml)

call set_calendar_type(calendar_type)

!-----------------------------------------------------------------------------!
! Read 3D-VAR GTS observation file
! at the end, 'ob' contains all observations
! iunit hardcoded at 99 to match 3DVAR_OBSPROC:DA_Setup_Obs_Structures.inc
!-----------------------------------------------------------------------------!
iunit = 99
open(iunit, file=gts_file)

max_sound_input = max_sound
max_synop_input = max_synop
max_satob_input = max_satob
max_airep_input = max_airep
max_satem_input = max_satem
max_pilot_input = max_pilot
max_amdar_input = max_amdar
max_metar_input = max_metar
max_gpspw_input = max_gpspw
max_ships_input = max_ships
max_ssmt1_input = max_ssmt1
max_ssmt2_input = max_ssmt2
max_ssmi_input  = max_ssmi
max_tovs_input  = max_tovs
max_qscat_input = max_qscat
max_profl_input = max_profl
max_buoys_input = max_buoys
max_bogus_input = max_bogus
max_gpsref_input = max_gpsref

CALL DA_Setup_Obs_Structures

!-----------------------------------------------------------------------------!
! Create DART obs sequence file
!-----------------------------------------------------------------------------!

! Initialize the counters:

num_obs = 0
num_obs_last = 0

call static_init_obs_sequence()

num_copies = 1
num_qc = 1

! Initialize an obs_sequence structure
call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)

call set_copy_meta_data(seq, 1, 'WRFVAR 2.1 GTS observations')
call set_qc_meta_data(  seq, 1, 'WRFVAR OBSPROC Quality Control')

call init_obs(obs, num_copies, num_qc)

!  -----------------
!  LOOP OVER RECORDS
!  -----------------

101 format(A,3I8)

!--- process sound  ---
if (Use_SoundObs .eqv. .TRUE.) then
   do i = 1, ob%num_sound
      call insert_gts_obs(ob%sound(i), &
                          qc_threshold = gts_qc_threshold, &
                          u_wind_type      = RADIOSONDE_U_WIND_COMPONENT, &
                          v_wind_type      = RADIOSONDE_V_WIND_COMPONENT, &
                          temperature_type = RADIOSONDE_TEMPERATURE, &
                          dew_point_type   = RADIOSONDE_DEWPOINT, &
                          which_vert=VERTISPRESSURE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
   write(*,101) 'Processed  OBS_TYPE    VAR#   DART#  Total_DART#'  
   write(*,101) 'Processed      TEMP', ob%num_sound, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process profiler  ---
if (Use_ProflObs .eqv. .TRUE.) then
   do i = 1, ob%num_profl
      call insert_gts_obs(ob%profl(i), &
                          qc_threshold = gts_qc_threshold, &
                          u_wind_type      = PROFILER_U_WIND_COMPONENT, &
                          v_wind_type      = PROFILER_V_WIND_COMPONENT, &
                          which_vert=VERTISPRESSURE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_profl, ' PROFI,  total obs # ', num_obs
   write(*,101) 'Processed     PROFI', ob%num_profl, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!---  process satob  ---
if (Use_SatobObs .eqv. .TRUE.) then
   do i = 1, ob%num_satob, max(1,num_thin_satob)      ! thinning: take one every num_thin_satob
      call insert_gts_obs(ob%satob(i), &
                          qc_threshold = gts_qc_threshold, &
                          u_wind_type      = SAT_U_WIND_COMPONENT, &
                          v_wind_type      = SAT_V_WIND_COMPONENT, &
                          which_vert=VERTISPRESSURE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_satob, ' SATOB,  total obs # ', num_obs
   write(*,101) 'Processed     SATOB', ob%num_satob, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process qscat ---
if (Use_QscatObs .eqv. .TRUE.) then
   do i = 1, ob%num_qscat, max(1,num_thin_qscat)      ! thinning: take one every num_thin_qscat
      call insert_gts_obs(ob%qscat(i), &
                          qc_threshold = gts_qc_threshold, &
                          u_wind_type      = QKSWND_U_WIND_COMPONENT, &
                          v_wind_type      = QKSWND_V_WIND_COMPONENT, &
                          which_vert=VERTISSURFACE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_qscat, ' QSCAT,  total obs # ', num_obs
   write(*,101) 'Processed     QSCAT', ob%num_qscat, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process synop ---
if (Use_SynopObs .eqv. .TRUE.) then
   do i = 1, ob%num_synop
      call insert_gts_obs(ob%synop(i), &
                          qc_threshold = gts_qc_threshold, &
                          pressure_type    = SYNOP_SURFACE_PRESSURE, &
                          u_wind_type      = SYNOP_U_WIND_COMPONENT, &
                          v_wind_type      = SYNOP_V_WIND_COMPONENT, &
                          temperature_type = SYNOP_TEMPERATURE, &
                          dew_point_type   = SYNOP_DEWPOINT, &
                          which_vert=VERTISSURFACE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_synop, ' SYNOP,  total obs # ', num_obs
   write(*,101) 'Processed     SYNOP', ob%num_synop, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process metar ---
if (Use_MetarObs .eqv. .TRUE.) then
   do i = 1, ob%num_metar
      call insert_gts_obs(ob%metar(i), &
                          qc_threshold = gts_qc_threshold, &
                          pressure_type    = METAR_SURFACE_PRESSURE, &
                          u_wind_type      = METAR_U_10_METER_WIND, &
                          v_wind_type      = METAR_V_10_METER_WIND, &
                          temperature_type = METAR_TEMPERATURE_2_METER, &
                          dew_point_type   = METAR_DEWPOINT_2_METER, &
                          which_vert=VERTISSURFACE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_metar, ' METAR,  total obs # ', num_obs
   write(*,101) 'Processed     METAR', ob%num_metar, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process buoys ---
if (Use_BuoysObs .eqv. .TRUE.) then
   do i = 1, ob%num_buoys
      call insert_gts_obs(ob%buoys(i), &
                          qc_threshold = gts_qc_threshold, &
                          pressure_type    = BUOY_SURFACE_PRESSURE, &
                          u_wind_type      = BUOY_U_WIND_COMPONENT, &
                          v_wind_type      = BUOY_V_WIND_COMPONENT, &
                          temperature_type = BUOY_TEMPERATURE, &
                          dew_point_type   = BUOY_DEWPOINT, &
                          which_vert=VERTISSURFACE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_buoys, ' BUOY ,  total obs # ', num_obs
   write(*,101) 'Processed     BUOYS', ob%num_buoys, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process ships ---
if (Use_ShipsObs .eqv. .TRUE.) then
   do i = 1, ob%num_ships
      call insert_gts_obs(ob%ships(i), & 
                          qc_threshold = gts_qc_threshold, &
                          pressure_type    = SHIP_SURFACE_PRESSURE, &
                          u_wind_type      = SHIP_U_WIND_COMPONENT, &
                          v_wind_type      = SHIP_V_WIND_COMPONENT, &
                          temperature_type = SHIP_TEMPERATURE, &
                          dew_point_type   = SHIP_DEWPOINT, &
                          which_vert=VERTISSURFACE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_ships, ' SHIP ,  total obs # ', num_obs
   write(*,101) 'Processed     SHIPS', ob%num_ships, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process airep ---
if (Use_AirepObs .eqv. .TRUE.) then
   do i = 1, ob%num_airep
      call insert_gts_obs(ob%airep(i), &
                          qc_threshold = gts_qc_threshold, &
                          u_wind_type      = AIREP_U_WIND_COMPONENT, & 
                          v_wind_type      = AIREP_V_WIND_COMPONENT, &
                          temperature_type = AIREP_TEMPERATURE, &
                          dew_point_type   = AIREP_DEWPOINT, &
                          which_vert=VERTISPRESSURE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_airep, ' AIREP,  total obs # ', num_obs
   write(*,101) 'Processed     AIREP', ob%num_airep, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif 

!--- process amdar ---
if (Use_AmdarObs .eqv. .TRUE.) then
   do i = 1, ob%num_amdar
      call insert_gts_obs(ob%amdar(i), &
                          qc_threshold = gts_qc_threshold, &
                          u_wind_type      = AMDAR_U_WIND_COMPONENT, &         
                          v_wind_type      = AMDAR_V_WIND_COMPONENT, &
                          temperature_type = AMDAR_TEMPERATURE, &
                          dew_point_type   = AMDAR_DEWPOINT, &
                          which_vert=VERTISPRESSURE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_amdar, ' AMDAR,  total obs # ', num_obs
   write(*,101) 'Processed     AMDAR', ob%num_amdar, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process pilot ---
if (Use_PilotObs .eqv. .TRUE.) then
   do i = 1, ob%num_pilot
      call insert_gts_obs(ob%pilot(i), &
                          qc_threshold = gts_qc_threshold, &
                          u_wind_type      = PILOT_U_WIND_COMPONENT, &
                          v_wind_type      = PILOT_V_WIND_COMPONENT, &
                          temperature_type = PILOT_TEMPERATURE, &
                          dew_point_type   = PILOT_DEWPOINT, &
                          which_vert=VERTISPRESSURE, num_obs=num_obs, obs=obs, seq=seq)
   enddo                  
!  print*,'Processed ', ob%num_pilot, ' PILOT,  total obs # ', num_obs
   write(*,101) 'Processed     PILOT', ob%num_pilot, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process bogus ---
if (Use_BogusObs .eqv. .TRUE.) then
   do i = 1, ob%num_bogus
      call insert_gts_obs(ob%bogus(i), &   
                          qc_threshold = gts_qc_threshold, &
                          u_wind_type      = BOGUS_U_WIND_COMPONENT, &
                          v_wind_type      = BOGUS_V_WIND_COMPONENT, &
                          temperature_type = BOGUS_TEMPERATURE, &
                          dew_point_type   = BOGUS_DEWPOINT, &
                          which_vert=VERTISPRESSURE, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_pilot, ' BOGUS,  total obs # ', num_obs
   write(*,101) 'Processed     BOGUS', ob%num_bogus, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif

!--- process gpsref ---
if (Use_gpsrefObs .eqv. .TRUE.) then
   gpskey = 0
   do i = 1, ob%num_gpsref
      call insert_gts_obs(ob%gpsref(i), &
                          qc_threshold = gts_qc_threshold, &
                          gpsref_type = GPSRO_REFRACTIVITY, &
                          gpsref_key = gpskey, gpsref_form = 1, isgps=.true.,&
                          which_vert=VERTISHEIGHT, num_obs=num_obs, obs=obs, seq=seq)
   enddo
!  print*,'Processed ', ob%num_gpsref, 'GPSREF,  total obs # ', num_obs
   write(*,101) 'Processed    GPSREF', ob%num_gpsref, num_obs-num_obs_last, num_obs
   num_obs_last = num_obs
endif


!!--- process satem  ---
!if (Use_SatemObs .eqv. .TRUE.) then
!!   do i = 1, ob%num_satem
!      call insert_gts_obs(ob%satem(i), &
!                          qc_threshold = gts_qc_threshold, &
!                          thickness_type   = SATEM_THICKNESS, &
!                          which_vert=VERTISPRESSURE, num_obs=num_obs, obs=obs, seq=seq)
!   enddo
!   write(*,101) 'Processed  OBS_TYPE    VAR#   DART#  Total_DART#'
!   write(*,101) 'Processed      TEMP', ob%num_satem, num_obs-num_obs_last, num_obs
!   num_obs_last = num_obs
!endif


call close_file(iunit)

call write_obs_seq(seq, obs_seq_out_file_name)

call error_handler(E_MSG, 'gts_to_dart', 'FINISHED gts_to_dart.')
call error_handler(E_MSG, 'gts_to_dart', 'Finished successfully.',&
                   source,revision,revdate)
call finalize_utilities()

 
END PROGRAM gts_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program make_fake_obs

! Generate an observation sequence FINAL file to thoroughly test obs_diag.
!
! Lots of goals for this routine.
! * all observation 'levels' for all 'VERTISxxxx' possibilities
! * exercise all 'regions'
! * exercise U,V wind combinations
!
! I want to test the vertical binning, the geographic binning, etc.
!

use         types_mod, only : r8, MISSING_R8, metadatalength
use     utilities_mod, only : initialize_utilities, get_unit, error_handler, E_ERR, &
                              open_file, close_file, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read, nmlfileunit, &
                              do_nml_file, do_nml_term
use  time_manager_mod, only : time_type, set_calendar_type, set_date, get_time
use      location_mod, only : VERTISUNDEF,    VERTISSURFACE, VERTISLEVEL, &
                              VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT
use obs_utilities_mod, only : add_obs_to_seq, create_3d_obs
use  obs_sequence_mod, only : obs_sequence_type, obs_type, init_obs, get_num_obs, &
                              static_init_obs_sequence, init_obs_sequence, &
                              set_copy_meta_data, set_qc_meta_data, write_obs_seq, &
                              set_obs_values, set_qc
use      obs_kind_mod, only : PRESSURE,                   &
                              ARGO_SALINITY,              &
                              SPECIFIC_HUMIDITY,          &
                              DROPSONDE_SURFACE_PRESSURE, &
                              ACARS_U_WIND_COMPONENT,     &
                              ACARS_V_WIND_COMPONENT,     &
                              EFGWORO

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2, string3
integer, parameter :: stringlength = 32

integer, parameter :: index_obs              = 1 ! 'observations'
integer, parameter :: index_truth            = 2 ! 'truth'
integer, parameter :: index_prior_ens_mean   = 3 ! 'prior ensemble mean'
integer, parameter :: index_poste_ens_mean   = 4 ! 'posterior ensemble mean'
integer, parameter :: index_prior_ens_spread = 5 ! 'prior ensemble spread'
integer, parameter :: index_poste_ens_spread = 6 ! 'posterior ensemble spread'
integer, parameter :: index_qc_in            = 1
integer, parameter :: index_qc_dart          = 2

integer, parameter :: MaxLevels       = 12
integer, parameter :: MaxRegions      = 4
integer, parameter :: MaxVertCoords   = 7
integer, parameter :: MaxObsPerRegion = 5
integer, parameter :: max_obs    = MaxRegions * MaxObsPerRegion * MaxLevels * MaxVertCoords
integer, parameter :: ens_size   = 10
integer, parameter :: num_copies = 6 + ens_size*2
integer, parameter :: num_qc     = 2

! 1) PRESSURE                     MODEL LEVEL    obs_def_gps_mod.f90
! 2) ARGO_SALINITY                HEIGHT         obs_def_ocean_mod.f90
! 3) SPECIFIC_HUMIDITY            SCALE HEIGHT   obs_def_gps_mod.f90
! 4) DROPSONDE_SURFACE_PRESSURE   SURFACE        obs_def_reanalysis_bufr_mod.f90
! 5) ACARS_U_WIND_COMPONENT       PRESSURE       obs_def_reanalysis_bufr_mod.f90
! 6) ACARS_V_WIND_COMPONENT       PRESSURE       obs_def_reanalysis_bufr_mod.f90
! 7) EFGWORO                      UNDEFINED      obs_def_GWD_mod.f90

integer, parameter, dimension(MaxVertCoords) :: vertlevels = &
      (/ VERTISLEVEL       , &
         VERTISHEIGHT      , &
         VERTISSCALEHEIGHT , &
         VERTISSURFACE     , &
         VERTISPRESSURE    , &
         VERTISPRESSURE    , &
         VERTISUNDEF     /)

integer, parameter, dimension(MaxVertCoords) :: obs_type_list = &
      (/ PRESSURE                   , &
         ARGO_SALINITY              , &
         SPECIFIC_HUMIDITY          , &
         DROPSONDE_SURFACE_PRESSURE , &
         ACARS_U_WIND_COMPONENT     , &
         ACARS_V_WIND_COMPONENT     , &
         EFGWORO /)

integer, dimension(MaxVertCoords) :: numlevels
integer  :: Nplevels=0, Nhlevels=0, Nmlevels=0

!-------------------------------------------------------------------------------
! SOME of the obs_diag namelist is required to determine observation locations.
!-------------------------------------------------------------------------------

real(r8), dimension(MaxLevels)   :: plevel       = MISSING_R8 ! pressure level (hPa)
real(r8), dimension(MaxLevels)   :: hlevel       = MISSING_R8 ! height (meters)
real(r8), dimension(MaxLevels)   :: mlevel       = MISSING_R8 ! model level (integer index)
real(r8), dimension(MaxLevels+1) :: plevel_edges = MISSING_R8 ! pressure level (hPa)
real(r8), dimension(MaxLevels+1) :: hlevel_edges = MISSING_R8 ! height (meters)
real(r8), dimension(MaxLevels+1) :: mlevel_edges = MISSING_R8 ! model levels (nondimensional)

real(r8), dimension(MaxRegions)  :: lonlim1 = (/   0.0,   0.0,   0.0, 235.0 /)
real(r8), dimension(MaxRegions)  :: lonlim2 = (/ 360.0, 360.0, 360.0, 295.0 /)
real(r8), dimension(MaxRegions)  :: latlim1 = (/  20.0, -80.0, -20.0,  25.0 /)
real(r8), dimension(MaxRegions)  :: latlim2 = (/  80.0, -20.0,  20.0,  55.0 /)

integer,  dimension(6)  :: first_bin_center = (/ 2003, 1, 1, 0, 0, 0 /)
integer,  dimension(6)  :: last_bin_center  = (/ 2003, 1, 2, 0, 0, 0 /)
integer,  dimension(6)  :: bin_separation   = (/    0, 0, 0, 6, 0, 0 /)
integer,  dimension(6)  :: bin_width        = (/    0, 0, 0, 6, 0, 0 /)
integer,  dimension(6)  :: time_to_skip     = (/    0, 0, 1, 0, 0, 0 /)

! Intentionally using only a subset of variables and then NOT calling
! the DART namelist error-checking routine  check_namelist_read()

namelist /obs_diag_nml/ plevel, hlevel, mlevel, &
                        plevel_edges, hlevel_edges, mlevel_edges, &
                        first_bin_center, last_bin_center, bin_separation, &
                        bin_width, time_to_skip

!-------------------------------------------------------------------------------
! variables associated with the observation
!-------------------------------------------------------------------------------

type(obs_type)          :: prev_obs, obs, next_obs
type(obs_sequence_type) :: seq
type(time_type)         :: prev_time, time_obs

integer  :: iunit, io
integer  :: iobs, iregion, ilevel, verttype
integer  :: year, month, day, hour, minute, second, obsec, obday
real(r8) :: lat, lon, vert, oberr, obsval, obsqc

logical  :: first_obs = .true.

!===============================================================================
! Get the ball rolling
!===============================================================================

call initialize_utilities('make_fake_obs')
call set_calendar_type('Gregorian')

call find_namelist_in_file("input.nml", "obs_diag_nml", iunit)
read(iunit, nml = obs_diag_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_diag_nml")

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=obs_diag_nml)
if (do_nml_term()) write(    *      , nml=obs_diag_nml)

call SetPressureLevels(plevel, plevel_edges, Nplevels)
call SetHeightLevels(  hlevel, hlevel_edges, Nhlevels)
call SetModelLevels(   mlevel, mlevel_edges, Nmlevels)

write(*,*)'pressure levels     = ',plevel(      1:Nplevels)
write(*,*)'pressure interfaces = ',plevel_edges(1:Nplevels+1)
write(*,*)'height   levels     = ',hlevel(      1:Nhlevels)
write(*,*)'height   interfaces = ',hlevel_edges(1:Nhlevels+1)
write(*,*)'model    levels     = ',mlevel(      1:Nmlevels)
write(*,*)'model    interfaces = ',mlevel_edges(1:Nmlevels+1)

! order:      LEVEL HEIGHT SCALEHEIGHT SURFACE PRESSURE PRESSURE UNDEF
numlevels = (/ Nmlevels, Nhlevels, Nhlevels, 1, Nplevels, Nplevels, 1 /)

! only 1 observation time for now.
year   = first_bin_center(1)
month  = first_bin_center(2)
day    = first_bin_center(3)
hour   = first_bin_center(4)
minute = first_bin_center(5)
second = first_bin_center(6)

time_obs = set_date(year, month, day, hour, minute, second)
call get_time(time_obs, obsec, obday)

call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
call init_obs(next_obs, num_copies, num_qc)
first_obs = .true.
prev_time = time_obs

! create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(seq, num_copies, num_qc, max_obs)
call generate_copy_meta_data(ens_size, seq)
call generate_qc_meta_data( seq )

! order:      LEVEL HEIGHT SCALEHEIGHT SURFACE PRESSURE PRESSURE UNDEF
coords  : do verttype = 1,MaxVertCoords
levels  : do ilevel   = 1,numlevels(verttype)
!regions : do iregion  = 1,MaxRegions
regions : do iregion  = 1,1
obsloop : do iobs     = 1,MaxObsPerRegion

   call CompleteObservation(vertlevels(verttype), obs_type_list(verttype), ilevel, iregion, iobs)

enddo obsloop
enddo regions
enddo levels
enddo coords

if ( get_num_obs(seq) > 0 ) then
   print *, 'writing obs_seq, obs_count = ', get_num_obs(seq)
   call write_obs_seq(seq, 'make_fake_obs.final')
endif

call finalize_utilities()

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------



Subroutine generate_copy_meta_data(ens_size, seq)

! This mimics what is done in filter:filter_generate_copy_meta_data()

integer,                 intent(in)    :: ens_size
type(obs_sequence_type), intent(inout) :: seq

integer :: i, copycount
character(len=metadatalength) :: prior_meta_data, posterior_meta_data

call set_copy_meta_data(seq, 1, 'observations')
call set_copy_meta_data(seq, 2, 'truth')
call set_copy_meta_data(seq, 3, 'prior ensemble mean')
call set_copy_meta_data(seq, 4, 'posterior ensemble mean')
call set_copy_meta_data(seq, 5, 'prior ensemble spread')
call set_copy_meta_data(seq, 6, 'posterior ensemble spread')

copycount = 6

! Set up obs ensemble members as requested
do i = 1, ens_size
   write(    prior_meta_data, '(a, 1x, i6)') 'prior ensemble member', i
   write(posterior_meta_data, '(a, 1x, i6)') 'posterior ensemble member', i
   copycount = copycount + 1
   call set_copy_meta_data(seq, copycount, prior_meta_data)
   copycount = copycount + 1
   call set_copy_meta_data(seq, copycount, posterior_meta_data)
end do

end Subroutine generate_copy_meta_data



Subroutine CompleteObservation( vertType, obsType, levelindx, regionindx, obsindx )

integer, intent(in) :: vertType    ! [VERTISUNDEF ... VERTISPRESSURE]
integer, intent(in) :: obsType     ! [PRESSURE ... EFGWORO]
integer, intent(in) :: levelindx   ! [1 .... N]
integer, intent(in) :: regionindx  ! [1 .... 4]
integer, intent(in) :: obsindx     ! [1 .... MaxObsPerRegion]

integer :: i, copycount

real(r8), dimension(1) :: obs_value, true_value, dartqc
real(r8), dimension(1) :: prior, prior_mean, prior_spread
real(r8), dimension(1) :: poste, poste_mean, poste_spread

call FillLocation (regionindx, obsindx, lon, lat )

! just put the ob right in the middle ... or on one edge ... or ...
if (    vertType == VERTISPRESSURE    ) then
!       vert     = plevel(levelindx)
        vert     = plevel_edges(levelindx)
elseif (vertType == VERTISLEVEL       ) then
!       vert     = mlevel(levelindx)
        vert     = mlevel_edges(levelindx)
elseif (vertType == VERTISHEIGHT      ) then
!       vert     = hlevel(levelindx)
        vert     = hlevel_edges(levelindx)
elseif (vertType == VERTISSCALEHEIGHT ) then
!       vert     = hlevel(levelindx)
        vert     = hlevel_edges(levelindx)
else
        vert     = 3.0_r8
endif

obsval = 1.0_r8
oberr  = 2.0_r8
obsqc  = 0.0_r8

call create_3d_obs(lat, lon, vert, vertType, obsval, &
             obsType,  oberr, obday, obsec, obsqc, obs)

obs_value(1)    = 10.0
true_value(1)   =  9.0
prior_mean(1)   =  8.0
poste_mean(1)   =  8.5
prior_spread(1) =  1.0
poste_spread(1) =  0.5

call set_obs_values(obs,    obs_value, 1)
call set_obs_values(obs,   true_value, 2)
call set_obs_values(obs,   prior_mean, 3)
call set_obs_values(obs,   poste_mean, 4)
call set_obs_values(obs, prior_spread, 5)
call set_obs_values(obs, poste_spread, 6)

copycount = 6

! Set up obs ensemble members as requested
do i = 1, ens_size
   prior(1) = real(i,r8)
   poste(1) = real(i,r8)
   copycount = copycount + 1
   call set_obs_values(obs, prior, copycount)
   copycount = copycount + 1
   call set_obs_values(obs, poste, copycount)
end do

! Lets make at least one ob have a bad DART value.
if (obsindx == 5) then
   dartqc(1) =  7.0_r8
else
   dartqc(1) =  0.0_r8
endif
call set_qc(obs, dartqc, index_qc_dart)

call add_obs_to_seq(seq, obs, time_obs, prev_obs, prev_time, first_obs)

end Subroutine CompleteObservation



Subroutine generate_qc_meta_data( seq )

type(obs_sequence_type), intent(inout) :: seq

call set_qc_meta_data(seq, 1, 'Quality Control')
call set_qc_meta_data(seq, 2, 'DART quality control')

end Subroutine generate_qc_meta_data




Subroutine FillLocation (regionindx, obsindx, lon, lat )
integer,  intent(in)  :: regionindx, obsindx
real(r8), intent(out) :: lon, lat

real(r8) :: dlat, dlon

dlat = (latlim2(regionindx) - latlim1(regionindx))/10.0_r8
dlon = (lonlim2(regionindx) - lonlim1(regionindx))/10.0_r8

if (    obsindx == 1)  then ! Upper left corner of region
   lat  = latlim2(regionindx) - dlat
   lon  = lonlim1(regionindx) + dlon
elseif (obsindx == 2)  then ! Upper right corner of region
   lat  = latlim2(regionindx) - dlat
   lon  = lonlim2(regionindx) - dlon
elseif (obsindx == 2)  then ! lower right corner of region
   lat  = latlim1(regionindx) + dlat
   lon  = lonlim2(regionindx) - dlon
elseif (obsindx == 4)  then ! lower left corner of region
   lat  = latlim1(regionindx) + dlat
   lon  = lonlim1(regionindx) + dlon
else                        ! dead center
   lat  = (latlim1(regionindx) + latlim2(regionindx)) / 2.0_r8
   lon  = (lonlim1(regionindx) + lonlim2(regionindx)) / 2.0_r8
endif
end Subroutine FillLocation



Subroutine SetPressureLevels(layerMiddles, layerEdges, nLayers)

! This is a utility routine ... presumes user input is minimalist
! and correct. Decreasing, no duplicate values ... etc.

real(r8), dimension(:), intent(inout) :: layerMiddles
real(r8), dimension(:), intent(inout) :: layerEdges
integer,                intent(  out) :: nLayers

integer :: i, n

! Count number of layer midpoints specified in namelist
n = 0
LevelLoop : do i = 1,MaxLevels
   if (layerMiddles(i) == MISSING_R8 ) exit LevelLoop
   n = n + 1
enddo LevelLoop

if ((n == 0) .or. (n == MaxLevels) ) then
   write(string1,*)n,' is too many plevel values - max is ',MaxLevels,' (MaxLevels)'
   call error_handler(E_ERR,'SetPressureLevels', string1,source,revision,revdate)
endif

! Set the layer edges ...
   layerEdges(  1) = layerMiddles(1)   - (layerMiddles(2) - layerMiddles( 1 ))/2.0_r8
   layerEdges(n+1) = layerMiddles(n)   + (layerMiddles(n) - layerMiddles(n-1))/2.0_r8
do i = 2,n
   layerEdges( i ) = layerMiddles(i-1) - (layerMiddles(i) - layerMiddles(i-1))/2.0_r8
enddo

layerEdges( 1 ) = min(1025.0_r8, layerEdges( 1 ))
layerEdges(n+1) = max(   0.0_r8, layerEdges(n+1))

nLayers = n

end Subroutine SetPressureLevels



Subroutine SetHeightLevels(layerMiddles, layerEdges, nLayers)

! This is a utility routine ... presumes user input is minimalist
! and correct. Ascending, no duplicate values ... etc.

real(r8), dimension(:), intent(inout) :: layerMiddles
real(r8), dimension(:), intent(inout) :: layerEdges
integer,                intent(  out) :: nLayers

integer :: i, n

! Count number of layer midpoints specified in namelist
n = 0
LevelLoop : do i = 1,MaxLevels
   if (layerMiddles(i) == MISSING_R8 ) exit LevelLoop
   n = n + 1
enddo LevelLoop

if ((n == 0) .or. (n == MaxLevels) ) then
   write(string1,*)n,' is too many hlevel values - max is ',MaxLevels,' (MaxLevels)'
   call error_handler(E_ERR,'SetHeightLevels', string1,source,revision,revdate)
endif

! Set the layer edges ...
   layerEdges(  1) = layerMiddles(1)   - (layerMiddles(2) - layerMiddles( 1 ))/2.0_r8
   layerEdges(n+1) = layerMiddles(n)   + (layerMiddles(n) - layerMiddles(n-1))/2.0_r8
do i = 2,n
   layerEdges( i ) = layerMiddles(i-1) + (layerMiddles(i) - layerMiddles(i-1))/2.0_r8
enddo

! FIXME ... what about depths ...
! layerEdges(1) = max(0.0_r8, layerEdges(1))

nLayers = n

end Subroutine SetHeightLevels



Subroutine SetModelLevels(layerMiddles, layerEdges, nLayers)

! This is a utility routine ... presumes user input is minimalist
! and correct. Ascending, no duplicate values ... etc.

real(r8), dimension(:), intent(inout) :: layerMiddles
real(r8), dimension(:), intent(inout) :: layerEdges
integer,                intent(  out) :: nLayers

integer :: i, n

! Count number of layer midpoints specified in namelist
n = 0
LevelLoop : do i = 1,MaxLevels
   if (layerMiddles(i) == MISSING_R8 ) exit LevelLoop
   n = n + 1
enddo LevelLoop

if ((n == 0) .or. (n == MaxLevels) ) then
   write(string1,*)n,' is too many mlevel values - max is ',MaxLevels,' (MaxLevels)'
   call error_handler(E_ERR,'SetModelLevels', string1,source,revision,revdate)
endif

! Set the layer edges ...
   layerEdges(  1) = layerMiddles(1)   - (layerMiddles(2) - layerMiddles( 1 ))/2.0_r8
   layerEdges(n+1) = layerMiddles(n)   + (layerMiddles(n) - layerMiddles(n-1))/2.0_r8
do i = 2,n
   layerEdges( i ) = layerMiddles(i-1) + (layerMiddles(i) - layerMiddles(i-1))/2.0_r8
enddo

! FIXME ... what about depths ...
layerEdges(1) = max(0.0_r8, layerEdges(1))

nLayers = n

end Subroutine SetModelLevels




end program make_fake_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

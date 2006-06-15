! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_def_radar_mod

! BEGIN DART PREPROCESS KIND LIST
! DOPPLER_RADIAL_VELOCITY, KIND_VELOCITY
! RADAR_REFLECTIVITY, KIND_RADAR_REFLECTIVITY
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_radar_mod, only : write_rad_vel, read_rad_vel, &
!                                 interactive_rad_vel, get_expected_rad_vel
!   use obs_def_radar_mod, only : get_expected_rad_ref
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(DOPPLER_RADIAL_VELOCITY)
!            call get_expected_rad_vel(state, location, obs_def%key, obs_val, istatus)
!         case(RADAR_REFLECTIVITY)
!            call get_expected_rad_ref(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(DOPPLER_RADIAL_VELOCITY)
!         call read_rad_vel(obs_def%key, ifile, fileformat)
!      case(RADAR_REFLECTIVITY)
!         continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(DOPPLER_RADIAL_VELOCITY)
!         call write_rad_vel(obs_def%key, ifile, fileformat)
!      case(RADAR_REFLECTIVITY)
!         continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(DOPPLER_RADIAL_VELOCITY)
!         call interactive_rad_vel(obs_def%key)
!      case(RADAR_REFLECTIVITY)
!         continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use        types_mod, only : r8, missing_r8, ps0, PI, gravity
use    utilities_mod, only : register_module, error_handler, E_ERR
use     location_mod, only : location_type, write_location, read_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                             KIND_TEMPERATURE, KIND_VERTICAL_VELOCITY, &
                             KIND_RAINWATER_MIXING_RATIO, KIND_DENSITY, &
                             KIND_GRAUPEL_MIXING_RATIO, KIND_SNOW_MIXING_RATIO

implicit none
private

public :: write_rad_vel, read_rad_vel, set_rad_vel, interactive_rad_vel, &
          get_expected_rad_vel, get_expected_rad_ref

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

logical, save :: module_initialized = .false.

! Storage for the special information required for observations of this type
real(r8), parameter :: dief = 0.224_r8
real(r8), parameter :: n0r = 8.0e6_r8, n0g = 4.0e4_r8, n0s = 3.0e6_r8
real(r8), parameter :: rho_r = 1000.0_r8, rho_g = 917.0_r8, rho_s = 100.0_r8
integer,  parameter :: max_rad_vel_obs = 300000
type(location_type) :: rad_loc(max_rad_vel_obs)
real(r8)            :: direction(3,max_rad_vel_obs)
integer :: keycount = 0 ! cumulative index into rad_loc,direction

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

subroutine write_rad_vel(key, ifile, fform)

integer,          intent(in)           :: key, ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      call write_location(ifile,  rad_loc(key), fileformat)
      call write_orientation(ifile, direction(:,key), fileformat)
      ! Write out the obs_def key for this observation
      write(ifile) key
   CASE DEFAULT
      ! Write the 5 character identifier for verbose formatted output
      write(ifile, 11)
11    format('platform')
      call write_location(ifile,  rad_loc(key), fileformat)
      call write_orientation(ifile, direction(:,key), fileformat)
      ! Write out the obs_def key for this observation
      write(ifile, *) key
END SELECT

end subroutine write_rad_vel



 subroutine read_rad_vel(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_rad_vel(key, ifile, fform)
!
! There is a whopping great list of metadata for each GPS obs.
! if you read multiple obs_sequence files, you need to keep
! appending into the (module) storage for 'rad_loc(:)' and 
! 'direction(:)' metadata streams. This is the difference
! between 'key' and 'keycount'

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=8)    :: header
character(len=32)   :: fileformat
real(r8)            :: orientation(3)
type(location_type) :: location

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      location = read_location(ifile, fileformat)
      orientation = read_orientation(ifile, fileformat)
      ! Read in the key for this particular observation
      read(ifile) key  ! basically, read and throw away
   CASE DEFAULT
      ! Read the character identifier for verbose formatted output
      read(ifile, FMT='(a8)') header
      if(header /= 'platform') then
         call error_handler(E_ERR,'read_rad_vel', &
              'Expected location header "platform" in input file', &
              source, revision, revdate)
      endif
      location = read_location(ifile, fileformat)
      orientation = read_orientation(ifile, fileformat)
      ! Read in the key for this particular observation
      read(ifile, *) key ! basically, read and throw away
END SELECT

keycount = keycount + 1    ! the total metadata key count from all sequences
key = keycount             ! copied to the output variable

if(key > max_rad_vel_obs) then
   write(*, *) 'key (',key,') exceeds max_rad_vel_obs (',max_rad_vel_obs,')'
   call error_handler(E_ERR,'read_rad_vel', &
              'Increase max_rad_vel_obs.', source, revision, revdate)
endif

rad_loc(key) = location
direction(:, key) = orientation

end subroutine read_rad_vel

!----------------------------------------------------------------------

subroutine set_rad_vel(key, location, orientation)

integer,             intent(in) :: key
real(r8),            intent(in) :: orientation(3)
type(location_type), intent(in) :: location

if ( .not. module_initialized ) call initialize_module

if(key > max_rad_vel_obs) then
   write(*, *) 'key (',key,') exceed max_rad_vel_obs (',max_rad_vel_obs,')'
   call error_handler(E_ERR,'read_rad_vel', &
              'Increase max_rad_vel_obs.', source, revision, revdate)
endif

rad_loc(key) = location
direction(:, key) = orientation

end subroutine set_rad_vel

!----------------------------------------------------------------------

subroutine interactive_rad_vel(key)

integer, intent(out) :: key

! Initializes the specialized part of a rad_vel observation
! Passes back up the key for this one

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(key >= max_rad_vel_obs) then
   call error_handler(E_ERR,'interactive_rad_vel', &
              'not ready to use', source, revision, revdate)
endif

! Increment the index
key = key + 1

end subroutine interactive_rad_vel

!----------------------------------------------------------------------

subroutine get_expected_rad_vel(state_vector, location, key, vr, istatus)

! Reference: Lin et al., 1983 (J. Climate Appl.Meteor., 1065-1092)
! Note that the reflectivity-weighted mean terminal velocities are used here.

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: vr
integer,             intent(out) :: istatus

real(r8), parameter :: a = 8.42e20_r8, b = 0.8_r8, c = 4.84e18_r8, d = 0.25_r8
real(r8), parameter :: CD = 0.6_r8
real(r8), parameter :: rhos0 = 1.0_r8
real(r8), parameter :: e = 4.0_r8*gravity*rho_g/(3.0_r8*CD), f = 0.5_r8
real(r8), parameter :: gam7b = 3376.92_r8, gam7d = 1155.38_r8, gam7f = 1871.25_r8
real(r8), parameter :: powr = (7.0_r8 + b)/4.0_r8
real(r8), parameter :: pows = (7.0_r8 + d)/4.0_r8
real(r8), parameter :: powg_dry = (7.0_r8 + f)/4.0_r8
real(r8), parameter :: powg_wet = 1.7875_r8

real(r8) :: u, v, w, qr, qg, qs, alpha, wt, rho, temp, precip, ref
real(r8) :: ar, as_wet, as_dry, ag_dry, ag_wet

if ( .not. module_initialized ) call initialize_module

call interpolate(state_vector, location, KIND_U_WIND_COMPONENT, u, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_V_WIND_COMPONENT, v, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_VERTICAL_VELOCITY, w, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_RAINWATER_MIXING_RATIO, qr, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_GRAUPEL_MIXING_RATIO, qg, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_SNOW_MIXING_RATIO, qs, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_DENSITY, rho, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_TEMPERATURE, temp, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif

precip = rho * (qr + qs + qg)
if(precip <= epsilon(precip)) then
   wt = 0.0_r8
else
   alpha=sqrt(rhos0/rho)
   ar = n0r*a*gam7b/(PI*rho_r*n0r)**powr
   wt = ar*((rho * qr)**powr)
   as_wet = n0s*c*gam7d/(PI*rho_s*n0s)**pows
   if ( temp < 273.15_r8 ) then
      as_dry = dief*((rho_s/rho_r)**2)*as_wet
      wt = alpha*(wt + as_dry*((rho * qs)**pows))
      ag_dry = 1.0e18_r8*dief*((rho_g/rho_r)**2)*n0g*gam7f/(PI*rho_g*n0g)**powg_dry
      wt = wt + sqrt(e/rho)*ag_dry*((rho * qg)**powg_dry)
   else
      wt = alpha*(wt + as_wet*((rho * qs)**pows))
      ag_wet = ((7.2e20_r8)**0.95_r8)*gam7f/(720.0_r8*(n0g**0.8375_r8))
      wt = wt + sqrt(e/rho)*ag_wet*(((rho * qg)/(PI*rho_g))**powg_wet)
   endif

   call get_reflectivity(qr, qg, qs, rho, temp, ref)
   wt = wt/ref

!!$   if(precip < epsilon(precip)) then
!!$      print*,'Terminal velocity = ',wt,qr,qs,qg,temp,10.0_r8 * log10(ref),epsilon(precip),precip
!!$   endif

endif

! direction(1) = sin(az)cos(elev)
! direction(2) = cos(az)cos(elev)
! direction(3) = sin(elev)
! az and elev are angles at the observation location.

vr = direction(1, key)*u + direction(2, key)*v + direction(3, key)*(w-wt)

end subroutine get_expected_rad_vel

!----------------------------------------------------------------------

subroutine get_expected_rad_ref(state_vector, location, ref, istatus)
!
! Computes "radar reflectivity" [ = 10 * log_10( radar reflectivity factor) ]
! in dBZ.  Reflectivities below ref_thresh are set to ref_thresh.
!
! ref_thresh hardwired.  Should be set in namelist.  [CS 17 Sep 04]

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: ref
integer,             intent(out) :: istatus

real(r8), parameter :: ref_thresh = 1.0_r8 ! Reflectivity below this value
                                           ! (corresponding to 0 dBZ) is set to this value.
                                           ! In future, would like to flag this occurence (through
                                           ! istatus??). Should be set in namelist.

real(r8) :: qr, qg, qs, rho, temp

if ( .not. module_initialized ) call initialize_module

call interpolate(state_vector, location, KIND_RAINWATER_MIXING_RATIO , qr , istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_GRAUPEL_MIXING_RATIO, qg, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_SNOW_MIXING_RATIO, qs, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_DENSITY, rho, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_TEMPERATURE, temp, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call get_reflectivity(qr, qg, qs, rho, temp, ref)

!!$ref = 10.0_r8 * log10(max(ref_thresh, ref))

end subroutine get_expected_rad_ref

!----------------------------------------------------------------------

subroutine get_reflectivity(qr, qg, qs, rho, temp, ref)
!
! Computes "radar reflectivity factor" in mm^6 m^-3
!
! References: Ferrier, 1994 (JAS, 249-280)
!             Smith et al., 1984 (JCAM, 1258-1260)
!             Smith at al., 1975 (JAM, 1156-1165)

real(r8), intent(in)  :: qr, qg, qs, rho, temp
real(r8), intent(out) :: ref

real(r8) :: precip

!!$According to Smith 1984 (JCAM, 1258-1260), there are two choices for
!!$the dielectric factor (dief), depending on how the snow particle sizes are specified.
!!$If melted raindrop diameters are used, then the factor is 0.224.  If
!!$equivalent ice sphere diameters are used, then the factor is 0.189.

!!$real(r8), parameter :: ar = 7.2e20_r8/(((PI*rho_r)**1.75_r8)*(n0r**0.75_r8))
!!$real(r8), parameter :: ag_dry = dief*((rho_g/rho_r)**2)*7.2e20_r8/ &
!!$                                (((PI*rho_g)**1.75_r8)*(n0g**0.75_r8))
!!$! This is appropriate for 10-cm radar.
!!$real(r8), parameter :: ag_wet = (7.2e20_r8/(((PI*rho_g)**1.75_r8)*(n0g**0.75_r8)))**0.95_r8
!!$
!!$real(r8), parameter :: as_wet = 7.2e20_r8/(((PI*rho_s)**1.75_r8)*(n0s**0.75_r8))
!!$real(r8), parameter :: as_dry = dief*((rho_s/rho_r)**2)*as_wet

!!$real(r8), parameter :: ar = 7.2e20_r8/((exp(log(PI*rho_r)*1.75_r8))*(exp(log(n0r)*0.75_r8)))
!!$real(r8), parameter :: ag_dry = dief*((rho_g/rho_r)**2)*7.2e20_r8/ &
!!$                                ((exp(log(PI*rho_g)*1.75_r8))*exp(log(n0g)*0.75_r8)))
!!$! This is appropriate for 10-cm radar.
!!$real(r8), parameter :: ag_wet = exp(log(7.2e20_r8/((exp(log(PI*rho_g)*1.75_r8))*(exp(log(n0g)*0.75_r8))))*0.95_r8)
!!$
!!$real(r8), parameter :: as_wet = 7.2e20_r8/((exp(log(PI*rho_s)*1.75_r8))*exp(log(n0s)*0.75_r8)))
!!$real(r8), parameter :: as_dry = dief*((rho_s/rho_r)**2)*as_wet

real(r8) :: ar, ag_dry, ag_wet, as_wet, as_dry

if ( .not. module_initialized ) call initialize_module

ref = 0.0_r8

! RAIN
precip = rho * qr
if ( precip > 0.0_r8 ) then
   ar = 7.2e20_r8/(((PI*rho_r)**1.75_r8)*(n0r**0.75_r8))
   ref = ref + ar * (precip**1.75_r8)
endif

! HAIL / GRAUPEL
precip = rho * qg
if ( precip > 0.0_r8 ) then
   if ( temp < 273.15_r8 ) then
      ag_dry = dief*((rho_g/rho_r)**2)*7.2e20_r8/ &
           (((PI*rho_g)**1.75_r8)*(n0g**0.75_r8))
      ref = ref + ag_dry * (precip**1.75_r8)
   else
      ! This is appropriate for 10-cm radar.
      ag_wet = (7.2e20_r8/(((PI*rho_g)**1.75_r8)*(n0g**0.75_r8)))**0.95_r8
      ref = ref + ag_wet * (precip**1.6625_r8)
   endif
endif

! SNOW
precip = rho * qs
if ( precip > 0.0_r8 ) then
   as_wet = 7.2e20_r8/(((PI*rho_s)**1.75_r8)*(n0s**0.75_r8))
   if ( temp < 273.15_r8 ) then
      as_dry = dief*((rho_s/rho_r)**2)*as_wet
      ref = ref + as_dry * (precip**1.75_r8)
   else
      ref = ref + as_wet * (precip**1.75_r8)
   endif
endif

end subroutine get_reflectivity

!----------------------------------------------------------------------

subroutine write_orientation(ifile, orientation, fform)

integer,                    intent(in) :: ifile
real(r8),                   intent(in) :: orientation(3)
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) orientation(1), orientation(2), orientation(3)
   CASE DEFAULT
      write(ifile, '(''dir3d'')' ) 
      write(ifile, *) orientation(1), orientation(2), orientation(3)
END SELECT

end subroutine write_orientation

!----------------------------------------------------------------------

function read_orientation(ifile, fform)

! Reads orientation from file that was written by write_orientation.
! See write_orientation for additional discussion.

integer,                    intent(in) :: ifile
real(r8)                               :: read_orientation(3)
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_orientation(1), read_orientation(2), read_orientation(3)
   CASE DEFAULT
      read(ifile, '(a5)' ) header

      if(header /= 'dir3d') then
         write(errstring,*)'Expected orientation header "dir3d" in input file, got ', header
         call error_handler(E_ERR, 'read_orientation', errstring, source, revision, revdate)
      endif
! Now read the orientation data value
      read(ifile, *) read_orientation(1), read_orientation(2), read_orientation(3)
END SELECT

end function read_orientation

!----------------------------------------------------------------------------

!!$subroutine simul_radar(max_num_obs, obs_num, num_copies, num_qc, &
!!$     obs_sequence, obs)
!!$
!!$integer,                 intent(in)    :: max_num_obs, num_copies, num_qc
!!$integer,                 intent(inout) :: obs_num
!!$type(obs_type),          intent(inout) :: obs
!!$type(obs_sequence_type), intent(inout) :: obs_sequence
!!$
!!$type(location_type)    :: rad_loc
!!$
!!$real(r8)               :: elev_clear(9)  = (/0.5_r8, 1.5_r8, 2.4_r8, 3.4_r8, &
!!$                                             4.3_r8, 6.0_r8, 9.9_r8, 14.6_r8, 19.5_r8/)
!!$real(r8)               :: elev_storm(14) = (/0.5_r8, 1.5_r8, 2.4_r8, 3.4_r8, &
!!$                                             4.3_r8, 5.3_r8, 6.2_r8, 7.5_r8, 8.7_r8, &
!!$                                             10.0_r8, 12.0_r8, 14.0_r8, 16.7_r8, 19.5_r8/)
!!$real(r8), allocatable  :: elev(:), azim(:), gate(:)
!!$real(r8)               :: faz, laz, daz, raz, fgate, lgate, dgate, rgate, var, elev_rad
!!$
!!$integer :: i, ilev, n_elev, iaz, n_az, igate, n_gate
!!$
!!$character(len=129) :: msgstring
!!$
!!$! Does interactive initialization of radar observations
!!$
!!$write(*, *)'Input radar location'
!!$call interactive_location(rad_loc)
!!$
!!$n_elev = 0
!!$
!!$do while(n_elev /= 9 .and. n_elev /= 14 .and. n_elev /= -1)
!!$   write(*, *)'Input 9 for 9 pre-defined elevations (clear mode)'
!!$   write(*, *)'Input 14 for 14 pre-defined elevations (storm mode)'
!!$   write(*, *)'Input -1 to choose number and elevation angles'
!!$   read(*, *) n_elev
!!$end do
!!$
!!$if(n_elev == -1) then
!!$   do while(n_elev < 1)
!!$      write(*, *)'Input number of elevations'
!!$      read(*, *) n_elev
!!$   end do
!!$   allocate(elev(n_elev))
!!$   do i = 1, n_elev
!!$      write(*, FMT='(a,i2)') 'Input elevation angle # ',i
!!$      read(*, *) elev(i)
!!$   end do
!!$else
!!$   allocate(elev(n_elev))
!!$   if(n_elev == 9) elev(:) = elev_clear(:)
!!$   if(n_elev == 14) elev(:) = elev_storm(:)
!!$endif
!!$
!!$write(*, *)'Input first azimuth angle (degree)'
!!$read(*, *) faz
!!$faz = DEG2RAD*faz
!!$
!!$write(*, *)'Input last azimuth angle (degree)'
!!$read(*, *) laz
!!$laz = DEG2RAD*laz
!!$
!!$write(*, *)'Input azimuth angle increment (degree)'
!!$read(*, *) daz
!!$daz = DEG2RAD*daz
!!$
!!$n_az = int((laz - faz)/daz) + 1
!!$
!!$allocate(azim(n_az))
!!$do iaz = 1, n_az
!!$   azim(iaz) = faz + (iaz-1)*daz
!!$enddo
!!$
!!$write(*, *)'Input closest gate (m)'
!!$read(*, *) fgate
!!$
!!$write(*, *)'Input farthest gate (m)'
!!$read(*, *) lgate
!!$
!!$write(*, *)'Input gate length (m)'
!!$read(*, *) dgate
!!$
!!$n_gate = int((lgate - fgate)/dgate) + 1
!!$
!!$allocate(gate(n_gate))
!!$do igate = 1, n_gate
!!$   gate(igate) = fgate + (igate-1)*dgate
!!$enddo
!!$
!!$write(*, *)'Input error variance for Doppler velocity'
!!$read(*, *) var
!!$
!!$do ilev = 1, n_elev
!!$
!!$   elev_rad = DEG2RAD*elev(ilev)
!!$
!!$   do iaz = 1, n_az
!!$
!!$      raz = azim(iaz)
!!$
!!$      do igate = 1, n_gate
!!$
!!$         rgate = gate(igate)
!!$
!!$         ! Radial velocity
!!$
!!$         if (obs_num <= max_num_obs) then
!!$
!!$            call set_radar_obs_def(rad_loc,rgate,raz,elev_rad,var,obs%def)
!!$
!!$            do i = 1, num_copies
!!$               write(*, *) 'Enter value ', i, 'for this observation'
!!$               read(*, *) obs%values(i)
!!$            end do
!!$
!!$            do i = 1, num_qc
!!$               write(*, *) 'Enter quality control value ', i, 'for this observation'
!!$               read(*, *) obs%qc(i)
!!$            end do
!!$
!!$            if(obs_num == 1) then
!!$               call insert_obs_in_seq(obs_sequence, obs)
!!$            else
!!$               call insert_obs_in_seq(obs_sequence, obs, &
!!$                    obs_sequence%obs(obs_num - 1))
!!$            endif
!!$
!!$            obs_num = obs_num + 1
!!$
!!$         else
!!$
!!$            write(msgstring,*) 'ran out of room, obs_num (',obs_num, &
!!$                 ') > max_num_obs (',max_num_obs,')'
!!$            call error_handler(E_ERR,'simul_radar',msgstring,source,revision,revdate)
!!$
!!$         endif
!!$
!!$         ! Reflectivity
!!$
!!$         if (obs_num <= max_num_obs) then
!!$
!!$            call set_radar_ref_obs_def(rad_loc,rgate,raz,elev_rad,var,obs%def)
!!$
!!$            do i = 1, num_copies
!!$               write(*, *) 'Enter value ', i, 'for this observation'
!!$               read(*, *) obs%values(i)
!!$            end do
!!$
!!$            do i = 1, num_qc
!!$               write(*, *) 'Enter quality control value ', i, 'for this observation'
!!$               read(*, *) obs%qc(i)
!!$            end do
!!$
!!$            if(obs_num == 1) then
!!$               call insert_obs_in_seq(obs_sequence, obs)
!!$            else
!!$               call insert_obs_in_seq(obs_sequence, obs, &
!!$                    obs_sequence%obs(obs_num - 1))
!!$            endif
!!$
!!$            obs_num = obs_num + 1
!!$
!!$         else
!!$
!!$            write(msgstring,*) 'ran out of room, obs_num (',obs_num, &
!!$                 ') > max_num_obs (',max_num_obs,')'
!!$            call error_handler(E_ERR,'simul_radar',msgstring,source,revision,revdate)
!!$
!!$         endif
!!$
!!$      end do
!!$   end do
!!$end do
!!$
!!$deallocate(gate, azim, elev)
!!$
!!$end subroutine simul_radar

!------------------------------------------------------------------------------

!!$subroutine set_radar_obs_def(rad_loc,rgate,raz,elev_rad,var,obs_def)
!!$
!!$! Allows creation of a radar radial velocity observation
!!$
!!$type(location_type),    intent(in)    :: rad_loc
!!$real(r8),               intent(in)    :: rgate,raz,elev_rad,var
!!$type(obs_def_type),     intent(inout) :: obs_def
!!$
!!$real(r8) :: h, spath, x, y, rad_lon, rad_lat, obs_lon, obs_lat, vloc
!!$real(r8) :: dir(3), elev_obs,ae
!!$integer  :: which_vert
!!$
!!$if ( .not. module_initialized ) call initialize_module
!!$
!!$! Set the observation kind
!!$obs_def%kind = KIND_VELOCITY
!!$
!!$vloc = query_location(rad_loc, 'vloc')
!!$which_vert = nint(query_location(rad_loc,'which_vert'))
!!$if(which_vert /= 3) then
!!$   call error_handler(E_ERR,'set_radar_obs_def', &
!!$        'Vertical coordinate of the radar location is not height', &
!!$        source, revision, revdate)
!!$endif
!!$
!!$! Doviak & Zrnic, 1993: Doppler radar and weather observations, eq. 2.28b-c
!!$
!!$ae = 4000.0_r8 * earth_radius / 3.0_r8
!!$
!!$h = sqrt( rgate*rgate + ae*ae + 2.0_r8 *rgate*ae*sin(elev_rad) ) - ae
!!$spath = ae * asin(rgate * cos(elev_rad) / (ae + h))
!!$
!!$x = spath*sin(raz)
!!$y = spath*cos(raz)
!!$
!!$rad_lon = query_location(rad_loc, 'lon')
!!$rad_lat = query_location(rad_loc, 'lat')
!!$
!!$ae = 1000.0_r8 * earth_radius
!!$
!!$obs_lat = y/ae + rad_lat
!!$obs_lon = x/(ae*cos(rad_lat + y/(2.0_r8*ae))) + rad_lon
!!$
!!$vloc = vloc + h
!!$
!!$obs_lon = obs_lon*RAD2DEG
!!$obs_lat = obs_lat*RAD2DEG
!!$
!!$obs_def%location = set_location(obs_lon, obs_lat, vloc, 3)
!!$
!!$! Set the time
!!$obs_def%time = set_time(0, 0)
!!$
!!$obs_def%error_variance = var
!!$
!!$! The elevation angle at the observation location is based on p. 23 of Battan's
!!$! radar meteorology book (equation 3.18), assuming an effective radius of the earth
!!$! of 4a/3.
!!$
!!$elev_obs = sqrt(1.5_r8 * h / (1000.0_r8 * earth_radius) + elev_rad*elev_rad)
!!$
!!$dir(1) = sin(raz)*cos(elev_obs)
!!$dir(2) = cos(raz)*cos(elev_obs)
!!$dir(3) = sin(elev_obs)
!!$
!!$!call set_platform_location(obs_def%platform, rad_loc)
!!$
!!$!call set_platform_orientation(obs_def%platform, dir)
!!$
!!$!call set_obs_def_platform(obs_def, obs_def%platform)
!!$
!!$end subroutine set_radar_obs_def
!!$
!!$!---------------------------------------------------------------------------
!!$
!!$subroutine set_radar_ref_obs_def(rad_loc,rgate,raz,elev_rad,var,obs_def)
!!$
!!$! Allows creation of a radar reflectivity observation
!!$
!!$type(location_type),    intent(in)    :: rad_loc
!!$real(r8),               intent(in)    :: rgate,raz,elev_rad,var
!!$type(obs_def_type),     intent(inout) :: obs_def
!!$
!!$real(r8) :: h, spath, x, y, rad_lon, rad_lat, obs_lon, obs_lat, vloc
!!$real(r8) :: ae
!!$
!!$if ( .not. module_initialized ) call initialize_module
!!$
!!$! Set the observation kind
!!$call set_obs_def_kind(obs_def,KIND_RADAR_REFLECTIVITY)
!!$
!!$vloc = query_location(rad_loc, 'vloc')
!!$
!!$! Doviak & Zrnic, 1993: Doppler radar and weather observations, eq. 2.28b-c
!!$
!!$ae = 4000.0_r8 * earth_radius / 3.0_r8
!!$
!!$h = sqrt( rgate*rgate + ae*ae + 2.0_r8 *rgate*ae*sin(elev_rad) ) - ae
!!$spath = ae * asin(rgate * cos(elev_rad) / (ae + h))
!!$
!!$x = spath*sin(raz)
!!$y = spath*cos(raz)
!!$
!!$rad_lon = query_location(rad_loc, 'lon')
!!$rad_lat = query_location(rad_loc, 'lat')
!!$
!!$ae = 1000.0_r8 * earth_radius
!!$
!!$obs_lat = y/ae + rad_lat
!!$obs_lon = x/(ae*cos(rad_lat + y/(2.0_r8*ae))) + rad_lon
!!$
!!$vloc = vloc + h
!!$
!!$obs_lon = obs_lon*RAD2DEG
!!$obs_lat = obs_lat*RAD2DEG
!!$
!!$obs_def%location = set_location(obs_lon, obs_lat, vloc, 3)
!!$
!!$! Set the time
!!$obs_def%time = set_time(0, 0)
!!$
!!$obs_def%error_variance = var
!!$
!!$end subroutine set_radar_ref_obs_def
!!$
!!$!----------------------------------------------------------------------------

end module obs_def_radar_mod

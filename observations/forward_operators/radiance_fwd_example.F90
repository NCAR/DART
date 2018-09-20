! example call to the forward operator for radiances.
! needs rttov libs, include files.
! needs to be preprocessed.

! Description:
!> @file
!!   Example calling direct model clear-sky simulation.
!
!> @brief
!!   Example calling direct model clear-sky simulation.
!!
!! @details
!!   This example is most easily run using the run_example_fwd.sh
!!   script (see the user guide).
!!
!!   This program requires the following files:
!!     the file containing input profiles (e.g. prof.dat)
!!     the RTTOV optical depth coefficient file
!!
!!   The output is written to a file called example_fwd_output.dat
!!
!!   This program demonstrates the most simple kind of RTTOV simulation.
!!
!!   You may wish to base your own code on this example in which case you
!!   should edit in particular the sections of code marked as follows:
!!       !================================
!!       !======Read =====start===========
!!            code to be modified
!!       !======Read ===== end ===========
!!       !================================
!!
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
PROGRAM example_fwd

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=11)  :: NameOfRoutine = 'example_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: prof_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  REAL(KIND=jprb)    :: trans_out(10)
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch
  INTEGER(KIND=jpim) :: np, nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV are as follows:
  !   1. Specify required RTTOV options
  !   2. Read coefficients
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof array with the channels/profiles to simulate
  !   5. Read input profile(s)
  !   6. Set up surface emissivity and/or reflectance
  !   7. Call rttov_direct and store results
  !   8. Deallocate all structures and arrays

  ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

!FIXME: things we need to set.
!  coef_filename comes from the rttov distribution
!  the profile comes from our forward operators
!  for now, number of profiles is 1
!  number of profile levels equals the number of model levels
!  no solar
!  pick 1 channel (need advice on this) - will be in obs in metadata
!  channel number is part of metadata
!  number of threads is 1 for now

  WRITE(0,*) 'enter path of coefficient file'
  READ(*,*) coef_filename
  WRITE(0,*) 'enter path of file containing profile data'
  READ(*,*) prof_filename
  WRITE(0,*) 'enter number of profiles'
  READ(*,*) nprof
  WRITE(0,*) 'enter number of profile levels'
  READ(*,*) nlevels
  WRITE(0,*) 'turn on solar simulations? (0=no, 1=yes)'
  READ(*,*) dosolar
  WRITE(0,*) 'enter number of channels to simulate per profile'
  READ(*,*) nchannels
  ALLOCATE(channel_list(nchannels))
  WRITE(0,*) 'enter space-separated channel list'
  READ(*,*,iostat=ios) channel_list(:)
  WRITE(0,*) 'enter number of threads to use'
  READ(*,*) nthreads


  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  IF (dosolar == 1) THEN
    opts % rt_ir % addsolar = .TRUE.           ! Include solar radiation
  ELSE
    opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
  ENDIF
  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method    ! FIXME = what are options?
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects
  opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects

  opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
  opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
  opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
  opts % rt_ir % co_data             = .FALSE. !
  opts % rt_ir % so2_data            = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !

  opts % config % verbose            = .TRUE.  ! Enable printing of warnings

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance, &
        init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchannels
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = channel_list(jch)
    ENDDO
  ENDDO


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============

! FIXME - call interpolate here for 

  OPEN(iup, file=TRIM(prof_filename), status='old', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening profile file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_skipcommentline(iup, errorstatus)

!FIXME - units for moisture, apparently ppmv or kg/kg
  ! Read gas units for profiles
  READ(iup,*) profiles(1) % gas_units
  profiles(:) % gas_units = profiles(1) % gas_units
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

    ! Read pressure (hPa), temp (K), WV, O3 (gas units ppmv or kg/kg - as read above)
    READ(iup,*) profiles(iprof) % p(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % t(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % q(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    ! Ozone profile is commented out in input profile data
!     READ(iup,*) profiles(iprof) % o3(:)
!     CALL rttov_skipcommentline(iup, errorstatus)

!FIXME - some models have surface t/q/p which are 2m
! but typically winds are 10m, and wind fetch isn't available.

    ! 2 meter air variables
    READ(iup,*) profiles(iprof) % s2m % t, &
                profiles(iprof) % s2m % q, &
                profiles(iprof) % s2m % p, &
                profiles(iprof) % s2m % u, &
                profiles(iprof) % s2m % v, &
                profiles(iprof) % s2m % wfetc
    CALL rttov_skipcommentline(iup, errorstatus)

!FIXME - might not have soil/water surface temp,
! don't have salinity in atm models, nor fastem for
! microwave emissivity models.
    ! Skin variables
    READ(iup,*) profiles(iprof) % skin % t,        &
                profiles(iprof) % skin % salinity, & ! Salinity only applies to FASTEM over sea
                profiles(iprof) % skin % fastem      ! FASTEM only applies to MW instruments
    CALL rttov_skipcommentline(iup, errorstatus)

!FIXME - don't have, in general
    ! Surface type and water type
    READ(iup,*) profiles(iprof) % skin % surftype, &
                profiles(iprof) % skin % watertype
    CALL rttov_skipcommentline(iup, errorstatus)

!FIXME - ok, these we understand.  verify elevation is in meters (km?)
    ! Elevation, latitude and longitude
    READ(iup,*) profiles(iprof) % elevation, &
                profiles(iprof) % latitude,  &
                profiles(iprof) % longitude
    CALL rttov_skipcommentline(iup, errorstatus)

!FIXME - these will be part of metadata.  note we generally
! have used azimuth/elevation instead of azimuth/zenith.
! elevation = 90 - zenith, and zenith = 90 - elevation.
    ! Satellite and solar angles
    READ(iup,*) profiles(iprof) % zenangle,    &
                profiles(iprof) % azangle,     &
                profiles(iprof) % sunzenangle, &
                profiles(iprof) % sunazangle
    CALL rttov_skipcommentline(iup, errorstatus)

!FIXME - make sure cfraction is 0 here
    ! Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
    READ(iup,*) profiles(iprof) % ctp, &
                profiles(iprof) % cfraction
    CALL rttov_skipcommentline(iup, errorstatus)

  ENDDO
  CLOSE(iup)

  !========== Read profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  ! In this example we have no values for input reflectances
  reflectance(:) % refl_in = 0._jprb

  ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
  ! (all channels in this case)
  calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)

  ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
  reflectance(:) % refl_cloud_top = 0._jprb


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_direct(                &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance) ! inout input/output BRDFs per channel
  ELSE
    CALL rttov_parallel_direct(     &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance,&! inout input/output BRDFs per channel
            nthreads    = nthreads)    ! in    number of threads to use
  ENDIF

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF


  !=====================================================
  !============== Output results == start ==============

  ! Open output file where results are written
  OPEN(ioout, file='output_'//NameOfRoutine//'.dat', status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' Instrument ', inst_name(coefs % coef % id_inst)
  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' '
  CALL rttov_print_opts(opts, lu=ioout)

  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    nprint = 1 + INT((nchannels-1)/10)
    WRITE(ioout,*)' '
    WRITE(ioout,*)' Profile ', iprof

    CALL rttov_print_profile(profiles(iprof), lu=ioout)

    WRITE(ioout,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
    WRITE(ioout,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
    WRITE(ioout,222) (radiance % bt(j), j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SATELLITE REFLECTANCES (BRF):'
      WRITE(ioout,444) (radiance % refl(j), j = 1+joff, nchannels+joff)
    ENDIF
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED RADIANCES (mW/m2/sr/cm-1):'
    WRITE(ioout,222) (radiance % total(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED OVERCAST RADIANCES:'
    WRITE(ioout,222) (radiance % cloudy(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE TO SPACE TRANSMITTANCE:'
    WRITE(ioout,4444) (transmission % tau_total(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE EMISSIVITIES:'
    WRITE(ioout,444) (emissivity(j) % emis_out, j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SURFACE BRDF:'
      WRITE(ioout,444) (reflectance(j) % refl_out, j = 1+joff, nchannels+joff)
    ENDIF

    IF (nchannels <= 20) THEN
      DO np = 1, nprint
          WRITE(ioout,*)' '
          WRITE(ioout,*)'Level to space transmittances for channels'
          WRITE(ioout,1115) (chanprof(j+joff) % chan, &
                    j = 1+(np-1)*10, MIN(np*10, nchannels))
          DO ilev = 1, nlevels
            DO j = 1 + (np-1)*10, MIN(np*10, nchannels)
              ! Select transmittance based on channel type (VIS/NIR or IR)
              IF (coefs % coef % ss_val_chn(chanprof(j+joff) % chan) == 2) THEN
                trans_out(j - (np-1)*10) = transmission % tausun_levels_path1(ilev,j+joff)
              ELSE
                trans_out(j - (np-1)*10) = transmission % tau_levels(ilev,j+joff)
              ENDIF
            ENDDO
            WRITE(ioout,4445) ilev, trans_out(1:j-1-(np-1)*10)
          ENDDO
          WRITE(ioout,1115) (chanprof(j+joff) % chan, &
                    j = 1+(np-1)*10, MIN(np*10, nchannels))
      ENDDO
    ENDIF
  ENDDO

  ! Close output file
  CLOSE(ioout, iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error closing the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  !============== Output results == end ==============
  !=====================================================


!FIXME - in the forward operator code we won't be deallocating
! this structure.

  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
1115 FORMAT(3X,10I8)
222  FORMAT(1X,10F8.2)
444  FORMAT(1X,10F8.3)
4444 FORMAT(1X,10F8.4)
4445 FORMAT(1X,I2,10F8.4)
777  FORMAT(/,A,A9,I3)

END PROGRAM example_fwd

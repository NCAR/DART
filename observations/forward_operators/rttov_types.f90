! Description:
!> @file
!!   Defines all derived types for RTTOV
!
!> @brief
!!   Defines all derived types for RTTOV
!!
!! @details
!!   This contains types that users will make use of in their code
!!   as well as types that RTTOV only uses internally.
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
MODULE rttov_types

  USE rttov_const, ONLY : &
      fastem_sp,          &
      ncldtyp,            &
      interp_rochon,      &
      ir_scatt_chou,      &
      vis_scatt_dom

  USE parkind1, ONLY : jpim, jprb, jplm

  IMPLICIT NONE

  ! ---------------------------------------------------------------------------
  ! User-level RTTOV structures
  ! ---------------------------------------------------------------------------

  !> Specify channels and profiles to simulate, declare as array of size
  !! nchanprof which is the total number of channels to simulate
  TYPE rttov_chanprof
    INTEGER(jpim) :: chan              !< Channel index
    INTEGER(jpim) :: prof              !< Profile index
  END TYPE

  !> Input/output surface emissivities, declare as array of size nchanprof
  TYPE rttov_emissivity
    REAL(jprb)    :: emis_in           !< Input emissivity (0-1)
    REAL(jprb)    :: emis_out          !< Output emissivity, value used by RTTOV (0-1)
  END TYPE

  !> Input/output surface BRDFs, declare as array of size nchanprof,
  !! can also be used to specify cloud top BRDF for simple cloud scheme in VIS/NIR channels
  TYPE rttov_reflectance
    REAL(jprb)    :: refl_in           !< Input BRDF (>=0)
    REAL(jprb)    :: refl_out          !< Output BRDF, value used by RTTOV (>=0)
    REAL(jprb)    :: refl_cloud_top    !< Optional, cloud top BRDF for simple cloud
  END TYPE

  !> Surface skin variables
  TYPE rttov_skin
    INTEGER(jpim) :: surftype          !< Surface type: 0=land, 1=sea, 2=sea-ice
    INTEGER(jpim) :: watertype         !< Water type: 0=fresh water, 1=ocean water
    REAL(jprb)    :: t                 !< Radiative skin temperature (K)
    REAL(jprb)    :: salinity          !< Practical ocean salinity unit (\%o) - FASTEM-4/5/6 only
    REAL(jprb)    :: foam_fraction     !< ocean foam fraction (0-1; only used if supply_foam_fraction is true)
    REAL(jprb)    :: snow_fraction     !< Snow coverage fraction for IR emissivity atlas (0-1)
    REAL(jprb)    :: soil_moisture     !< Soil moisture (m^3/m^3) - not currently used
    REAL(jprb)    :: fastem(fastem_sp) !< FASTEM land/sea-ice surface parameters
  END TYPE rttov_skin

  !> Surface 2m variables
  TYPE rttov_s2m
    REAL(jprb) :: t                    !< Temperature (K)
    REAL(jprb) :: q                    !< Water vapour (ppmv or kg/kg)
    REAL(jprb) :: o                    !< Ozone (ppmv or kg/kg) - not currently used
    REAL(jprb) :: p                    !< Surface pressure (hPa)
    REAL(jprb) :: u                    !< U 10m wind component (m/s)
    REAL(jprb) :: v                    !< V 10m wind component (m/s)
    REAL(jprb) :: wfetc                !< Wind fetch (metres)
  END TYPE rttov_s2m

  !> Atmospheric profiles on model pressure levels
  TYPE rttov_profile
    CHARACTER(LEN=128) :: id             !< Optional profile ID string
    INTEGER(jpim) :: date(3)             !< Year, Month, Day
    INTEGER(jpim) :: time(3)             !< Hour, Minute, Second

    INTEGER(jpim) :: nlevels             !< Number of atmospheric levels
    INTEGER(jpim) :: nlayers             !< Number of atmospheric layers

    INTEGER(jpim) :: gas_units           !< Units of gas profiles (0 or less => ppmv over dry air,
                                         !! 1 => kg/kg over moist air, 2 => ppmv over moist air)

    REAL(jprb), POINTER :: p(:)          !< Pressure (hPa)
    REAL(jprb), POINTER :: t(:)          !< Temperature (K)
    REAL(jprb), POINTER :: q(:)          !< Water vapour (ppmv or kg/kg)
    REAL(jprb), POINTER :: o3(:)         !< O3 (ppmv or kg/kg)
    REAL(jprb), POINTER :: co2(:)        !< CO2 (ppmv or kg/kg)
    REAL(jprb), POINTER :: n2o(:)        !< N2O (ppmv or kg/kg)
    REAL(jprb), POINTER :: co(:)         !< CO (ppmv or kg/kg)
    REAL(jprb), POINTER :: ch4(:)        !< CH4 (ppmv or kg/kg)
    REAL(jprb), POINTER :: so2(:)        !< SO2 (ppmv or kg/kg)
    REAL(jprb), POINTER :: clw(:)        !< Cloud liquid water absorption only (kg/kg)

    LOGICAL(jplm)       :: mmr_cldaer    !< Cloud/aerosol units: False => num density cm-3 (aer), g.m-3 (cld);
                                         !!                      True => kg/kg (cld, aer)
    REAL(jprb), POINTER :: aerosols(:,:) !< Aerosol layer concentrations (units: see mmr_cldaer)
    REAL(jprb), POINTER :: cloud(:,:)    !< Cloud liquid/ice water layer concentrations (units: see mmr_cldaer)
    REAL(jprb), POINTER :: cfrac(:)      !< Layer cloud fraction (0-1)
    REAL(jprb), POINTER :: clwde(:)      !< Cloud liquid water particle effective diameter (microns)
    INTEGER(jpim)       :: clw_scheme    !< Select liquid water cloud scheme (1-2)
    REAL(jprb), POINTER :: icede(:)      !< Ice particle effective diameter (microns)
    INTEGER(jpim)       :: idg           !< Ice particle effective diameter parameterisation (1-4)
    INTEGER(jpim)       :: ice_scheme    !< Select ice cloud scheme (1-2)

    TYPE(rttov_skin) :: skin             !< Surface skin variables
    TYPE(rttov_s2m)  :: s2m              !< Surface 2m variables

    REAL(jprb) :: zenangle               !< Satellite zenith angle (degrees)
    REAL(jprb) :: azangle                !< Satellite azimuth angle (degrees)
    REAL(jprb) :: sunzenangle            !< Solar azimuth angle (degrees)
    REAL(jprb) :: sunazangle             !< Solar azimuth angle (degrees)
    REAL(jprb) :: elevation              !< Surface elevation (km)
    REAL(jprb) :: latitude               !< Latitude (degrees)
    REAL(jprb) :: longitude              !< Longitude (degrees)

    REAL(jprb) :: Be                     !< Earth magnetic field strength (Gauss)
    REAL(jprb) :: cosbk                  !< Cosine of the angle between the Earth magnetic
                                         !!   field and wave propagation direction

    REAL(jprb) :: ctp                    !< Black body (simple) cloud top pressure (hPa)
    REAL(jprb) :: cfraction              !< Black body (simple) cloud fraction (0-1)
  END TYPE rttov_profile

  !> Additional atmospheric cloud/hydrometeor profile input for RTTOV-SCATT
  TYPE rttov_profile_cloud
    INTEGER(jpim) :: nlevels             !< Number of atmospheric levels (same as in rttov_profile)
    LOGICAL(jplm) :: use_totalice        !< False => separate ice and snow; True => total ice
    LOGICAL(jplm) :: mmr_snowrain        !< Snow and rain input units are: False => kg/m2/s; True => kg/kg
    REAL(jprb)    :: cfrac               !< Average cloud fraction (only used if lusercfrac = TRUE, 0-1)

    REAL(jprb), POINTER :: ph(:)         !< nlevels+1 of half-level model pressures (hPa)
    REAL(jprb), POINTER :: cc(:)         !< nlevels of cloud cover (0-1)
    REAL(jprb), POINTER :: clw(:)        !< nlevels of cloud liquid water (kg/kg)
    REAL(jprb), POINTER :: ciw(:)        !< nlevels of cloud ice water (kg/kg)
    REAL(jprb), POINTER :: totalice(:)   !< nlevels of total ice (kg/kg)
    REAL(jprb), POINTER :: rain(:)       !< nlevels of rain (units: see mmr_snowrain)
    REAL(jprb), POINTER :: sp(:)         !< nlevels of solid precipitation (units: see mmr_snowrain)
  END TYPE rttov_profile_cloud

  !> @internal Data for interpolating phase functions
  TYPE rttov_phasefn_int
    REAL(jprb)             :: zminphadiff
    REAL(jprb),    POINTER :: cosphangle(:)
    INTEGER(jpim), POINTER :: iphangle(:)
  END TYPE rttov_phasefn_int

  !> Explicit optical parameters for IR scattering
  TYPE rttov_opt_param
    REAL(jprb),    POINTER :: abs(:,:)       !< Absorption coef (nlayers,nchannels) (km-1)
    REAL(jprb),    POINTER :: sca(:,:)       !< Scattering coef (nlayers,nchannels) (km-1)
    REAL(jprb),    POINTER :: bpr(:,:)       !< b parameter (nlayers,nchannels) (no units)
    INTEGER(jpim)          :: nmom           !< Number of Leg. coefs provided for each phase fn
    REAL(jprb),    POINTER :: legcoef(:,:,:) !< Phase fn Leg. coefs (1:nmom+1,nlayers,nchannels)
                                             !! - required for DOM scattering calculations
    REAL(jprb),    POINTER :: pha(:,:,:)     !< Phase function (nphangle,nlayers,nchannels), should be
                                             !! normalised such that integral over all angles is 4*pi
                                             !! - required for solar scattering calculations
    REAL(jprb),    POINTER :: phangle(:)     !< Angles over which phase fns defined (nphangle) (degrees)
                                             !! - required for solar scattering calculations

    ! The following are for RTTOV internal purposes only
    TYPE(rttov_phasefn_int) :: phasefn_int   !< Data for interpolating phase functions - for internal use only
  END TYPE rttov_opt_param

  !> Output transmittances
  TYPE rttov_transmission
    REAL(jprb), POINTER  :: tau_total(:)             !< Surface-satellite transmittance (channels)
    REAL(jprb), POINTER  :: tau_levels(:,:)          !< Level-satellite transmittance (levels,channels)
    REAL(jprb), POINTER  :: tausun_total_path2(:)    !< Sun-surface-satellite solar transmittance
    REAL(jprb), POINTER  :: tausun_levels_path2(:,:) !< Sun-level-satellite solar transmittance for each level
    REAL(jprb), POINTER  :: tausun_total_path1(:)    !< Surface-satellite solar transmittance
    REAL(jprb), POINTER  :: tausun_levels_path1(:,:) !< Level-satellite solar transmittance for each level
  END TYPE rttov_transmission

  !> Output radiances and corresponding brightness temperatures and reflectances (BRFs)
  !! Radiance unit: mW/m2/sr/cm-1; BT unit: K; BRFs are unitless.
  !! Array sizes are (nchannels) or (nlayers,nchannels)
  TYPE rttov_radiance
    REAL(jprb),    POINTER :: clear(:)       !< Clear sky radiance
    REAL(jprb),    POINTER :: total(:)       !< Cloudy radiance for given cloud
    REAL(jprb),    POINTER :: bt_clear(:)    !< Brightness temp equivalent to clear radiance
    REAL(jprb),    POINTER :: bt(:)          !< Brightness temp equivalent to total radiance
    REAL(jprb),    POINTER :: refl_clear(:)  !< Reflectance calculated from clear radiance
    REAL(jprb),    POINTER :: refl(:)        !< Reflectance calculated from total radiance
    REAL(jprb),    POINTER :: overcast(:,:)  !< Overcast radiance for opaque cloud at level bounding
                                             !!   bottom of each layer
    REAL(jprb),    POINTER :: cloudy(:)      !< 100% cloudy radiance for given cloud (simple cloud scheme)
                                             !!   or same as total (addclouds/addaerosl true)
    INTEGER(jpim), POINTER :: quality(:)     !< Flag to indicate whether they may be accuracy issues with
                                             !!   simulated radiances
    LOGICAL(jplm)          :: plane_parallel !< Flag to indicate if strict plane-parallel was enforced
                                             !!   e.g. by user-specified option or for DOM
  END TYPE rttov_radiance

  !> Secondary radiances optionally calculated in direct model only for clear-sky with no solar contribution
  !! Radiance unit: mW/m2/sr/cm-1. Array sizes are (nchannels) or (nlayers,nchannels)
  TYPE rttov_radiance2
    REAL(jprb), POINTER  :: upclear(:)     !< Clear sky upwelling radiance without reflection term
    REAL(jprb), POINTER  :: dnclear(:)     !< Clear sky downwelling radiance
    REAL(jprb), POINTER  :: refldnclear(:) !< Reflected clear sky downwelling radiance
    REAL(jprb), POINTER  :: up(:,:)        !< Sum( B * dT ) above cloud upwelling radiance from each layer
    REAL(jprb), POINTER  :: down(:,:)      !< Sum( B / T**2 dT ) above cloud downwelling radiance from
                                           !!   each layer
    REAL(jprb), POINTER  :: surf(:,:)      !< Radiance at surface emitted from a black cloud
  END TYPE rttov_radiance2

  !> Output variables from RTTOV-SCATT enabling emissivity retrievals
  !! All arrays are size(chanprof)
  TYPE rttov_scatt_emis_retrieval_type
    REAL(jprb), POINTER  :: cfrac(:)    !< RTTOV-SCATT effective cloud fraction (Tallsky = cfrac * Tcld + (1-cfrac) * Tclr
    REAL(jprb), POINTER  :: bsfc(:)     !< Surface black-body radiance (i.e. Planck(Tsfc))
    REAL(jprb), POINTER  :: tau_cld(:)  !< Along-path transmittance, surface to space (cloudy column)
    REAL(jprb), POINTER  :: up_cld(:)   !< TOA upwelling radiance or TB from atmosphere (not inc. surface emission or reflection)
    REAL(jprb), POINTER  :: down_cld(:) !< SFC downwelling radiance or TB (inc. cosmic term)
    REAL(jprb), POINTER  :: tau_clr(:)  !< Along-path transmittance, surface to space (clear column)
    REAL(jprb), POINTER  :: up_clr(:)   !< TOA upwelling radiance or TB from atmosphere (not inc. surface emission or reflection)
    REAL(jprb), POINTER  :: down_clr(:) !< SFC downwelling radiance or TB (inc. cosmic term)
  END TYPE rttov_scatt_emis_retrieval_type

  !> Output PC scores and reconstructed radiances from PC-RTTOV
  !! Radiance unit: mW/m2/sr/cm-1; BT unit: K.
  TYPE rttov_pccomp
    REAL(jprb), POINTER  :: pcscores(:)    !< Principal component scores
    REAL(jprb), POINTER  :: total_pccomp(:)!< Radiances reconstructed using principal components
    REAL(jprb), POINTER  :: bt_pccomp(:)   !< Brightness temp equivalent to radiances
                                           !! reconstructed using principal components
  END TYPE rttov_pccomp


  !> Configuration options that apply to all flavours of RTTOV
  TYPE rttov_opts_config
    LOGICAL(jplm) :: apply_reg_limits = .FALSE.  !< Switch to restrict input profiles to coef training limits
    LOGICAL(jplm) :: verbose          = .TRUE.   !< Switch for verbose output
    LOGICAL(jplm) :: do_checkinput    = .TRUE.   !< Switch to apply internal profile checking
    LOGICAL(jplm) :: fix_hgpl         = .FALSE.  !< Switch to apply fix to match 2m p with elevation in geometry calculations
  END TYPE

  !> Options for PC-RTTOV
  TYPE rttov_opts_pc
    LOGICAL(jplm) :: addpc     = .FALSE.   !< Switch to enable PC-RTTOV
    INTEGER(jpim) :: ipcbnd    = -1        !< PC spectral band
    INTEGER(jpim) :: ipcreg    = -1        !< PC predictor channel set
    LOGICAL(jplm) :: addradrec = .FALSE.   !< Switch for calculation of reconstructed radiances
  END TYPE

  !> General radiative transfer options
  TYPE rttov_opts_rt_all
    LOGICAL(jplm) :: addrefrac      = .FALSE. !< Switch to enable atmospheric refraction
    LOGICAL(jplm) :: switchrad      = .FALSE. !< Switch for input units in AD/K models
    LOGICAL(jplm) :: use_q2m        = .TRUE.  !< Switch to enable use of 2m q variable
    LOGICAL(jplm) :: do_lambertian  = .FALSE. !< Switch for setting Lambertian reflection (IR and MW)
    LOGICAL(jplm) :: plane_parallel = .FALSE. !< Switch to ignore atmospheric curvature
    LOGICAL(jplm) :: dtau_test      = .TRUE.  !< Switch to apply dtau test in transmit/integrate calculations
  END TYPE

  !> VIS/IR-only radiative transfer options
  TYPE rttov_opts_rt_ir
    TYPE(rttov_opts_pc) :: pc                         !< PC-RTTOV options

    INTEGER(jpim) :: solar_sea_brdf_model = 1         !< Solar sea BRDF model (1-2)
    INTEGER(jpim) :: ir_sea_emis_model    = 2         !< IR sea emissivity model (1-2)
    LOGICAL(jplm) :: addsolar             = .FALSE.   !< Switch to enable solar simulations
    LOGICAL(jplm) :: do_nlte_correction   = .FALSE.   !< Switch to enable NLTE bias correction
    LOGICAL(jplm) :: addaerosl            = .FALSE.   !< Switch to enable IR aerosol calculations
    LOGICAL(jplm) :: addclouds            = .FALSE.   !< Switch to enable IR cloudy calculations
    LOGICAL(jplm) :: user_aer_opt_param   = .FALSE.   !< Switch to supply aerosol optical properties explicitly per channel
    LOGICAL(jplm) :: user_cld_opt_param   = .FALSE.   !< Switch to supply cloud optical properties explicitly per channel
    REAL(jprb)    :: cldstr_threshold     = -1.0_jprb !< Ignore cloud streams with weights lower than this
    LOGICAL(jplm) :: cldstr_simple        = .FALSE.   !< Switch for simplified cloud stream option (USE WITH CAUTION)

    INTEGER(jpim) :: ir_scatt_model  = ir_scatt_chou  !< IR scattering model to use
    INTEGER(jpim) :: vis_scatt_model = vis_scatt_dom  !< VIS/NIR scattering model to use

    INTEGER(jpim) :: dom_nstreams        = 8          !< Number of DOM streams, must be even and not less than 2
    REAL(jprb)    :: dom_accuracy        = 0._jprb    !< Convergence criterion for termination of DOM azimuthal loop
    REAL(jprb)    :: dom_opdep_threshold = 0._jprb    !< DOM ignores levels below this optical depth:
                                                      !!   10. is reasonable, not applied if <= 0

    LOGICAL(jplm) :: ozone_data = .FALSE.             !< Switch to enable input of O3 profile
    LOGICAL(jplm) :: co2_data   = .FALSE.             !< Switch to enable input of CO2 profile
    LOGICAL(jplm) :: n2o_data   = .FALSE.             !< Switch to enable input of N2O profile
    LOGICAL(jplm) :: co_data    = .FALSE.             !< Switch to enable input of CO profile
    LOGICAL(jplm) :: ch4_data   = .FALSE.             !< Switch to enable input of CH4 profile
    LOGICAL(jplm) :: so2_data   = .FALSE.             !< Switch to enable input of SO2 profile
  END TYPE

  !> MW-only radiative transfer options
  TYPE rttov_opts_rt_mw
    INTEGER(jpim) :: fastem_version        = 6        !< FASTEM version (0-6); 0 => TESSEM2
    LOGICAL(jplm) :: clw_data              = .FALSE.  !< Switch to enable input of cloud liquid water profile
    INTEGER(jpim) :: clw_scheme            = 1        !< MW CLW scheme: 1 => Liebe, 2 => Rosenkranz, 3 => TKC
    LOGICAL(jplm) :: clw_calc_on_coef_lev  = .TRUE.   !< Apply MW CLW calculations on coef/user levels (true/false resp.)
    LOGICAL(jplm) :: supply_foam_fraction  = .FALSE.  !< Supply a foam fraction to FASTEM
    LOGICAL(jplm) :: apply_band_correction = .TRUE.   !< Apply band-correction for Planck radiance and BT calculations
  END TYPE

  !> Options for internal vertical interpolation and vertical grid setup
  TYPE rttov_opts_interp
    LOGICAL(jplm) :: addinterp        = .FALSE.       !< Switch to enable RTTOV interpolator
    INTEGER(jpim) :: interp_mode      = interp_rochon !< Interpolation mode (1-5, see user guide)
    LOGICAL(jplm) :: lgradp           = .FALSE.       !< Switch to make pressure an active variable in TL/AD/K models
    LOGICAL(jplm) :: spacetop         = .TRUE.        !< Switch to assume space boundary at top-most input pressure level
    LOGICAL(jplm) :: reg_limit_extrap = .FALSE.       !< Switch to extrapolate input profiles using regression limits
  END TYPE

  !> HTFRTC options structure
  TYPE rttov_opts_htfrtc
    LOGICAL(jplm) :: htfrtc      = .FALSE. !< Switch to use htfrtc
    INTEGER(jpim) :: n_pc_in     = -1      !< Number of principal components to be used
    LOGICAL(jplm) :: reconstruct = .FALSE. !< Switch to select reconstructed radiances
  END TYPE rttov_opts_htfrtc

  !> Developer-only options structure
  TYPE rttov_opts_dev
    LOGICAL(jplm) :: no_opt_param_tladk    !< Ignore cld/aer_opt_param_tl/ad/k even if passed as arguments
  END TYPE rttov_opts_dev

  !> RTTOV options structure
  TYPE rttov_options
    TYPE(rttov_opts_config)  :: config          !< General configuration options
    TYPE(rttov_opts_rt_all)  :: rt_all          !< General RT options
    TYPE(rttov_opts_rt_ir)   :: rt_ir           !< VIS/IR RT options
    TYPE(rttov_opts_rt_mw)   :: rt_mw           !< MW RT options
    TYPE(rttov_opts_interp)  :: interpolation   !< Interpolation options
    TYPE(rttov_opts_htfrtc)  :: htfrtc_opts     !< HTFRTC options
    TYPE(rttov_opts_dev)     :: dev             !< Developer-only options
  END TYPE

  !> RTTOV-SCATT options structure: RTTOV-SCATT deliberately does
  !! not give user control over certain core RT options.
  TYPE rttov_options_scatt
    TYPE(rttov_opts_config) :: config                       !< General configuration options
    LOGICAL(jplm)  :: lusercfrac            = .FALSE.       !< Switch to enable user-specified effective cloud fraction
    LOGICAL(jplm)  :: apply_band_correction = .TRUE.        !< Apply band-correction for Planck radiance and BT calculations
    LOGICAL(jplm)  :: use_q2m               = .TRUE.        !< Switch to enable use of 2m q variable
    LOGICAL(jplm)  :: dtau_test             = .TRUE.        !< Switch to apply dtau test in transmit/integrate calculations
    INTEGER(jpim)  :: fastem_version        = 6             !< FASTEM version (1-6)
    LOGICAL(jplm)  :: supply_foam_fraction  = .FALSE.       !< Supply a foam fraction to FASTEM
    INTEGER(jpim)  :: interp_mode           = interp_rochon !< Interpolation mode (1-5, see user guide)
    LOGICAL(jplm)  :: reg_limit_extrap      = .FALSE.       !< Switch to extrapolate input profiles using regression limits
    LOGICAL(jplm)  :: lgradp                = .FALSE.       !< Switch to make pressure an active variable in TL/AD/K models
    LOGICAL(jplm)  :: lradiance             = .FALSE.       !< Switch for radiative transfer with radiances rather than TBs
  END TYPE

  ! ---------------------------------------------------------------------------
  ! Internal RTTOV structures
  ! ---------------------------------------------------------------------------

  !> @internal Satellite geometry
  TYPE rttov_geometry
    REAL(jprb) :: sinzen
    REAL(jprb) :: sinzen_sq
    REAL(jprb) :: coszen
    REAL(jprb) :: coszen_sq
    REAL(jprb) :: seczen
    REAL(jprb) :: seczen_sq
    REAL(jprb) :: sinview
    REAL(jprb) :: sinview_sq
    REAL(jprb) :: cosview_sq
    REAL(jprb) :: normzen
    REAL(jprb) :: coszen_sun
    REAL(jprb) :: sinzen_sun
    REAL(jprb) :: sinlat
    REAL(jprb) :: coslat
  END TYPE rttov_geometry

  !> @internal The actual predictor arrays
  TYPE rttov_path_pred
    REAL(jprb), POINTER :: mixedgas(:,:)   ! (nmixed,  nlayers)
    REAL(jprb), POINTER :: watervapour(:,:)! (nwater,  nlayers)
    REAL(jprb), POINTER :: ozone(:,:)      ! (nozone,  nlayers)
    REAL(jprb), POINTER :: wvcont(:,:)     ! (nwvcont, nlayers)
    REAL(jprb), POINTER :: co2(:,:)        ! (nco2,    nlayers)
    REAL(jprb), POINTER :: n2o(:,:)        ! (nn2o,    nlayers)
    REAL(jprb), POINTER :: co(:,:)         ! (nco,     nlayers)
    REAL(jprb), POINTER :: ch4(:,:)        ! (nch4,    nlayers)
    REAL(jprb), POINTER :: so2(:,:)        ! (nso2,    nlayers)
    REAL(jprb), POINTER :: pmc(:,:,:)      ! pressure modulated cell (npmc,nlevels,nchannels)
  END TYPE rttov_path_pred

  !> @internal Predictors
  TYPE rttov_predictors
    ! the nxxxx could be set to 0 to indicate the abscence
    ! of the predictor, in that case there is no need to
    ! allocate the corresponding predictor
    INTEGER(jpim) :: nlevels   ! number of levels for predictors (all same)
    INTEGER(jpim) :: nmixed    ! number of variables for Mixed Gases
    INTEGER(jpim) :: nwater    ! number of variables for Water Vapour
    INTEGER(jpim) :: nozone    ! number of variables for Ozone
    INTEGER(jpim) :: nwvcont   ! number of variables for WV Continuum
    INTEGER(jpim) :: nco2      ! number of variables for CO2
    INTEGER(jpim) :: nn2o      ! number of variables for N2O
    INTEGER(jpim) :: nco       ! number of variables for CO
    INTEGER(jpim) :: nch4      ! number of variables for CH4
    INTEGER(jpim) :: nso2      ! number of variables for SO2
    INTEGER(jpim) :: npmc      ! number of variables for pressure modulated cell correction

    TYPE(rttov_path_pred), POINTER :: path1(:)  ! Predictors for surface-satellite path (always required)
    TYPE(rttov_path_pred), POINTER :: path2(:)  ! Predictors for sun-surface-satellite path (only required for solar)
  END TYPE rttov_predictors

  !> @internal Actual storage for gas fast coefficients 
  TYPE rttov_fast_coef_gas
    REAL(jprb), POINTER :: coef(:,:) => NULL()
  END TYPE rttov_fast_coef_gas

  !> @internal Fast coefficients
  !  Separate structure to allow for non-PW and PW if necessary
  TYPE rttov_fast_coef
    ! SIZE is fmv_gas, i.e. one element per gas.
    TYPE(rttov_fast_coef_gas), POINTER  :: gasarray(:) => NULL()

    ! SHAPE is (ncoef, levels): these point to coefs within gasarray for ease of use
    REAL(jprb), POINTER :: mixedgas(:,:)    => NULL()
    REAL(jprb), POINTER :: watervapour(:,:) => NULL()
    REAL(jprb), POINTER :: ozone(:,:)       => NULL()
    REAL(jprb), POINTER :: wvcont(:,:)      => NULL()
    REAL(jprb), POINTER :: co2(:,:)         => NULL()
    REAL(jprb), POINTER :: n2o(:,:)         => NULL()
    REAL(jprb), POINTER :: co(:,:)          => NULL()
    REAL(jprb), POINTER :: ch4(:,:)         => NULL()
    REAL(jprb), POINTER :: so2(:,:)         => NULL()
  END TYPE rttov_fast_coef

  !> @internal RTTOV-11.2 style fast coefficient structure to maintain backwards
  !! compatibility for later versions
  TYPE rttov_fast_coef_hdf_io
    ! SHAPE is (levels, channels, ncoef)
    REAL(jprb), POINTER :: mixedgas(:,:,:)     => NULL()
    REAL(jprb), POINTER :: watervapour(:,:,:)  => NULL()
    REAL(jprb), POINTER :: ozone(:,:,:)        => NULL()
    REAL(jprb), POINTER :: wvcont(:,:,:)       => NULL()
    REAL(jprb), POINTER :: co2(:,:,:)          => NULL()
    REAL(jprb), POINTER :: n2o(:,:,:)          => NULL()
    REAL(jprb), POINTER :: co(:,:,:)           => NULL()
    REAL(jprb), POINTER :: ch4(:,:,:)          => NULL()
    REAL(jprb), POINTER :: so2(:,:,:)          => NULL()
  END TYPE rttov_fast_coef_hdf_io

  !> @internal NLTE coefficients
  TYPE rttov_nlte_coef
    REAL(jprb), POINTER :: coef(:,:,:,:) ! ncoef x nsat x nsol x nchan
    REAL(jprb), POINTER :: sol_zen_angle(:), sat_zen_angle(:)
    REAL(jprb), POINTER :: cos_sol(:), sec_sat(:)
    INTEGER(jpim)       :: ncoef, nsol, nsat, nchan
    INTEGER(jpim)       :: start_chan, end_chan
    REAL(jprb)          :: max_sat_angle
  END TYPE rttov_nlte_coef

  !> @internal Optical depth coefs
  TYPE rttov_coef
     ! Structure for the storage of RTTOV coefficients
     ! this may differ from what is stored in the coefficient files especially
     ! for the units (ie kg/kg to ppmv)
     ! Gases are separated in MxG WV O3
     ! Number of levels is the same for all gases (taken from MxG).
     !
    INTEGER(jpim)       :: id_platform         ! platform   (see user guide or rttov_const)
    INTEGER(jpim)       :: id_sat              ! satellite  (.....)
    INTEGER(jpim)       :: id_inst             ! instrument (.....)
    INTEGER(jpim)       :: id_sensor           ! 1 = Infrared, 2 = Microwave, 3 = High-resolution, 4 = Polarimeter
    INTEGER(jpim)       :: id_comp_lvl         ! RTTOV coefficient file version number
    INTEGER(jpim)       :: id_comp_pc          ! Principal component coefficient file version number
    INTEGER(jpim)       :: id_creation_date(3) ! YYYY MM DD
    CHARACTER (LEN=80)  :: id_creation         ! Creation comment
    CHARACTER (LEN=32)  :: id_common_name      ! usual name of the satellite
    CHARACTER (LEN=132) :: line_by_line(100)   ! readme LBL
    CHARACTER (LEN=132) :: readme_srf(100)     ! readme Spectral Response Function

       !FUNDAMENTAL_CONSTANTS section
    REAL(jprb)          :: fc_planck_c1        ! first radiation constant (mW/(m2*sr*cm-4))
    REAL(jprb)          :: fc_planck_c2        ! second radiation constant (cm*K)
    REAL(jprb)          :: fc_sat_height       ! satellite nominal altitude (km)

       !FAST_MODEL_VARIABLES section
    CHARACTER (LEN=32)          :: fmv_model_def  ! FMV definition (RTTOV6 OPTRAN RTTOV7)
    INTEGER(jpim)               :: fmv_model_ver  ! fast model version compatibility level
    INTEGER(jpim)               :: fmv_ori_nchn   ! number of channels in original file
    INTEGER(jpim)               :: fmv_chn        ! number of channels read from file into coef structure
    INTEGER(jpim)               :: fmv_gas        ! number of gases in file
    INTEGER(jpim), POINTER      :: fmv_gas_id(:)  ! gas id. number i gas_id list (fmv_gas)
    INTEGER(jpim), POINTER      :: fmv_gas_pos(:) ! respective position of each gas of gas_id list (ngases_max)
    INTEGER(jpim), POINTER      :: fmv_var(:)     ! number of variables/predictors by gas (fmv_gas)
    INTEGER(jpim), POINTER      :: fmv_coe(:)     ! number of coefficients by gas (fmv_gas)
    INTEGER(jpim), POINTER      :: fmv_lvl(:)     ! number of levels(pres/absorber) by gas (fmv_gas)

    INTEGER(jpim)               :: nmixed         ! number of variables/predictors for Mixed Gases
    INTEGER(jpim)               :: nwater         ! number of variables/predictors for Water Vapour
    INTEGER(jpim)               :: nozone         ! number of variables/predictors for Ozone
    INTEGER(jpim)               :: nwvcont        ! number of variables/predictors for WV continuum
    INTEGER(jpim)               :: nco2           ! number of variables/predictors for CO2
    INTEGER(jpim)               :: nn2o           ! number of variables/predictors for N2O
    INTEGER(jpim)               :: nco            ! number of variables/predictors for CO
    INTEGER(jpim)               :: nch4           ! number of variables/predictors for CH4
    INTEGER(jpim)               :: nso2           ! number of variables/predictors for SO2
    INTEGER(jpim)               :: nlevels        ! number of levels(pres/absorber) same for all gases
    INTEGER(jpim)               :: nlayers        ! number of layers(pres/absorber) nlevels-1
    LOGICAL(jplm)               :: inczeeman      ! Flag to include Zeeman effect for this sensor
    INTEGER(jpim)               :: ncmixed        ! number of coefficients for Mixed Gases
    INTEGER(jpim)               :: ncwater        ! number of coefficients for Water Vapour
    INTEGER(jpim)               :: ncozone        ! number of coefficients for Ozone
    INTEGER(jpim)               :: ncwvcont       ! number of coefficients for WV continuum
    INTEGER(jpim)               :: ncco2          ! number of coefficients for CO2
    INTEGER(jpim)               :: ncn2o          ! number of coefficients for N2O
    INTEGER(jpim)               :: ncco           ! number of coefficients for CO
    INTEGER(jpim)               :: ncch4          ! number of coefficients for CH4
    INTEGER(jpim)               :: ncso2          ! number of coefficients for SO2

       !FILTER_FUNCTIONS section  array size is fmv_chn
    LOGICAL(jplm)          :: ff_val_bc       ! are any band corrections to be applied?
    LOGICAL(jplm)          :: ff_val_gam      ! any gamma corrections?
    INTEGER(jpim), POINTER :: ff_ori_chn(:)   ! original chan number
    INTEGER(jpim), POINTER :: ff_val_chn(:)   ! validity of the channel (1=OK)
    REAL(jprb),    POINTER :: ff_cwn(:)       ! cental wave number (cm-1)
    REAL(jprb),    POINTER :: ff_bco(:)       ! band correction offset (K)
    REAL(jprb),    POINTER :: ff_bcs(:)       ! band correction slope (K/K)
    REAL(jprb),    POINTER :: ff_gam(:)       ! gamma factor transmittance correction

       !TRANSMITTANCE_TRESHOLD section  array size is fmv_chn
    INTEGER(jpim), POINTER :: tt_val_chn(:)
    REAL(jprb),    POINTER :: tt_a0(:)
    REAL(jprb),    POINTER :: tt_a1(:)

       !PLANCK_WEIGHTED section array size if fmv_chn
    INTEGER(jpim), POINTER :: pw_val_chn(:)   ! 0 => non-PW thermal coefs, 1 => PW thermal coefs

       !SOLAR_SPECTRUM section array size is fmv_chn
    INTEGER(jpim), POINTER :: ss_val_chn(:)
    REAL(jprb),    POINTER :: ss_solar_spectrum(:)  ! TOA solar irradiance (mW/m2/cm-1)
    REAL(jprb),    POINTER :: refl_visnir_ow(:)     ! Vis/NIR ocean water refl, populated by init_coef
    REAL(jprb),    POINTER :: refl_visnir_fw(:)     ! Vis/NIR fresh water refl, populated by init_coef

       !WATER_OPTICAL_CONSTANT section array size is fmv_chn
    COMPLEX(jprb), POINTER :: woc_waopc_ow(:)       ! Refractive index of ocean water
    COMPLEX(jprb), POINTER :: woc_waopc_fw(:)       ! Refractive index of fresh water

       !WAVE_SPECTRUM section array size is ws_nomega
       !Data used to compute the frequency spectrum of the JONSWAP
       !wave model surface wave.
    INTEGER(jpim)          :: ws_nomega
    REAL(jprb),    POINTER :: ws_npoint(:)
    REAL(jprb),    POINTER :: ws_k_omega(:)

       !FASTEM section
    INTEGER(jpim), POINTER :: fastem_polar(:)  ! polarisation of each channel
       ! 0 = 0.5 V+H
       ! 1 = 90 - incident angle
       ! 2 = incident angle
       ! 3 = vertical
       ! 4 = horizontal
       ! 5 = V+H
       ! Full stokes vector

       !SSIREM section     array size is fmv_chn
       ! ems =   ssirem_a0
       !       - ssirem_a1*(zen**ssirem_xzn1)
       !       - ssirem_a2*(zen**ssirem_xzn2)
       ! where zen is satellite zenith angle in degrees, divided by 60.
    INTEGER(jpim)          :: ssirem_ver      ! version number
    REAL(jprb),    POINTER :: ssirem_a0(:)    ! constant coef
    REAL(jprb),    POINTER :: ssirem_a1(:)    ! first coef
    REAL(jprb),    POINTER :: ssirem_a2(:)    ! second coef
    REAL(jprb),    POINTER :: ssirem_xzn1(:)  ! 1st exponent on zenith angle
    REAL(jprb),    POINTER :: ssirem_xzn2(:)  ! 2nd exponent on zenith angle

       ! IREMIS - IR sea surface emissivity model section
    INTEGER(jpim)         :: iremis_version   ! version number
    REAL(jprb)            :: iremis_angle0    ! reference zenith angle (degrees)
    REAL(jprb)            :: iremis_tskin0    ! reference Tskin (K)
    INTEGER(jpim)         :: iremis_ncoef     ! version number
    REAL(jprb),   POINTER :: iremis_coef(:,:) ! coefficients (iremis_ncoef,fmv_chn)

       !REFERENCE_PROFILE section  defined on mixed gases pressure levels
       ! gases are in the order of gas id codes
       ! unit for mr in coef file is ppmv over dry air
       ! unit for mr for optical depth calculations is ppmv over dry air
    REAL(jprb), POINTER      :: ref_prfl_p(:)     ! pressure  (hPa)       (levels)
    REAL(jprb), POINTER      :: ref_prfl_t(:,:)   ! temperature (K)       (levels, gases)
    REAL(jprb), POINTER      :: ref_prfl_mr(:,:)  ! mixing ratio (ppmv)   (levels, gases)
       !Background profiles to use when no optional input gas profile supplied
    REAL(jprb), POINTER      :: bkg_prfl_mr(:,:)  ! mixing ratio (ppmv)   (levels, gases)

       !PROFILE_ENVELOPE section
       ! gases are in the order of gas id codes
       ! unit for mr in coef file is ppmv over dry air
       ! unit for mr for optical depth calculations is ppmv over dry air
    REAL(jprb), POINTER      :: lim_prfl_p(:)       ! pressure  (hPa)       (levels)
    REAL(jprb), POINTER      :: env_prfl_tmax(:)    ! max temperature (K)   (levels)
    REAL(jprb), POINTER      :: env_prfl_tmin(:)    ! min temperature (K)   (levels)
    REAL(jprb), POINTER      :: env_prfl_gmax(:,:)  ! max mixing r (ppmv) (levels, gases)
    REAL(jprb), POINTER      :: env_prfl_gmin(:,:)  ! min mixing r (ppmv) (levels, gases)
       ! Profile limits calculated from envelope
    REAL(jprb), POINTER      :: lim_prfl_tmax(:)    ! max temperature (K)   (levels)
    REAL(jprb), POINTER      :: lim_prfl_tmin(:)    ! min temperature (K)   (levels)
    REAL(jprb), POINTER      :: lim_prfl_gmax(:,:)  ! max mixing r (ppmv) (levels, gases)
    REAL(jprb), POINTER      :: lim_prfl_gmin(:,:)  ! min mixing r (ppmv) (levels, gases)

       !FAST_COEFFICIENTS/SOLAR_FAST_COEFFICIENTS section
       ! For non-PW instruments, "solar" will point to "thermal" coefs
       ! For instruments with solar-affected PW channels, both thermal and solar
       ! structures will be populated from coef file
       ! Fast coefficients now grouped per channel rather than per gas (RTTOV12)
       ! so thermal and solar are now arrays with size nchannels.
    TYPE(rttov_fast_coef), POINTER :: thermal(:)          ! FAST_COEFFICIENTS
    TYPE(rttov_fast_coef), POINTER :: solar(:)            ! SOLAR_FAST_COEFFICIENTS when present
    LOGICAL(jplm)                  :: solarcoef           ! .TRUE. if solar fast coefs present in file: need
                                                          ! to know whether solar points to thermal or is
                                                          ! allocated separately.
    LOGICAL(jplm)                  :: nltecoef = .FALSE.  ! .TRUE. if nlte coefs present
    TYPE(rttov_nlte_coef), POINTER :: nlte_coef           ! nlte_coef

         !PRESSURE_MODULATED_CELL section
    LOGICAL(jplm)       :: pmc_shift = .FALSE. ! .TRUE. if pmc shift coefs present
    REAL(jprb)          :: pmc_lengthcell      ! cell length (cm)
    REAL(jprb), POINTER :: pmc_pnominal(:)     ! nominal cell pressure (hPa) - nchannels
    REAL(jprb)          :: pmc_tempcell        ! cell temperature (K)
    REAL(jprb)          :: pmc_betaplus1       ! co2 band-average: self-HW/air-HW
    INTEGER(jpim)       :: pmc_nlay            ! number of layers used
    INTEGER(jpim)       :: pmc_nvar            ! number of variables used
    REAL(jprb), POINTER :: pmc_coef(:,:,:)     ! pressure moodulated cell corrections - nlevels, nchannels, nvariables
    REAL(jprb), POINTER :: pmc_ppmc(:)         ! actual cell pressure (hPa) - nchannels

       ! Auxiliary variables
    REAL(jprb)               :: ratoe             ! ratio (H+R)/R  H=sat height, R=Earth radius
    INTEGER(jpim)            :: mwcldtop          ! Upper layer index for MW CLW absorber calcs
    REAL(jprb), POINTER      :: planck1(:)        ! C1 * Nu**3 (mW/(m2*sr*cm-1))
    REAL(jprb), POINTER      :: planck2(:)        ! C2 * Nu (K)
    REAL(jprb), POINTER      :: frequency_ghz(:)  ! frequency in GHz

       ! other predictor variables see Science and Validation report
    REAL(jprb), POINTER      :: dp(:)        ! interval between standard p levels (hPa)
    REAL(jprb), POINTER      :: dpp(:)       ! pressure based variable (hPa**2)
    REAL(jprb), POINTER      :: tstar(:)     ! layer temp (K)
    REAL(jprb), POINTER      :: to3star(:)   ! layer temp for O3 calculations (K)
    REAL(jprb), POINTER      :: wstar(:)     ! layer WV  (ppmv)
    REAL(jprb), POINTER      :: ostar(:)     ! layer O3  (ppmv)
    REAL(jprb), POINTER      :: co2star(:)   ! layer co2 (ppmv)
    REAL(jprb), POINTER      :: n2ostar(:)   ! layer n2o (ppmv)
    REAL(jprb), POINTER      :: costar(:)    ! layer co  (ppmv)
    REAL(jprb), POINTER      :: ch4star(:)   ! layer ch4 (ppmv)
    REAL(jprb), POINTER      :: so2star(:)   ! layer so2 (ppmv)

    INTEGER(jpim), POINTER   :: bounds(:,:,:,:)
  END TYPE rttov_coef

  !> @internal RTTOV-SCATT coefs
  TYPE rttov_scatt_coef
    ! Structure for the storage of RTTOV_SCATT coefficients
    INTEGER(jpim) :: nhydro ! Number of hydrometeors in computation
    INTEGER(jpim) :: mtype  ! Number of hydrometeors     in Mie tables
    INTEGER(jpim) :: mfreqm ! Number of frequencies      in Mie tables
    INTEGER(jpim) :: mtemp  ! Number of temperature bins in Mie tables
    INTEGER(jpim) :: mwc    ! Number of water bins       in Mie tables
    REAL(jprb)    :: offset_temp_rain       ! temperature offset in table for rain type (K)
    REAL(jprb)    :: offset_temp_sp         ! temperature offset in table for solid prec. type (K)
    REAL(jprb)    :: offset_temp_liq        ! temperature offset in table for cloud water type (K)
    REAL(jprb)    :: offset_temp_ice        ! temperature offset in table for cloud ice type (K)
    REAL(jprb)    :: offset_temp_totalice   ! temperature offset in table for total ice type (K)
    REAL(jprb)    :: offset_water           ! liquid/ice water offset in table
    REAL(jprb)    :: scale_water            ! log10(liquid/ice water) scaling factor in table
    REAL(jprb)    :: from_scale_water       ! 10**(1/scale_water)
    REAL(jprb)    :: conv_rain(2)           ! coefficients for rain unit conversion (mm.h-1 to g.m-3)
    REAL(jprb)    :: conv_sp  (2)           ! coefficients for solid prec. unit conversion (mm.h-1 to g.m-3)
    REAL(jprb)    :: conv_liq (2)           ! coefficients for cloud water conversion (not used)
    REAL(jprb)    :: conv_ice (2)           ! coefficients for cloud ice conversion   (not used)
    REAL(jprb)    :: conv_totalice (2)      ! coefficients for total ice conversion   (not used)
    REAL(jprb), POINTER :: mie_freq(:)      ! list of frequencies in Mie table
    REAL(jprb), POINTER :: ext(:,:,:,:)     ! extinction coefficent table
    REAL(jprb), POINTER :: ssa(:,:,:,:)     ! single scattering albedo table
    REAL(jprb), POINTER :: asp(:,:,:,:)     ! asymmetry parameter table
  END TYPE rttov_scatt_coef

  !> @internal Surface and cloud fraction
  TYPE rttov_profile_aux_s
    INTEGER(jpim) :: nearestlev_surf ! nearest model level above surface
    REAL(jprb)    :: pfraction_surf  ! pressure fraction of surface in model layer (hPa)
    INTEGER(jpim) :: nearestlev_ctp  ! nearest model level above cloud top
    REAL(jprb)    :: pfraction_ctp   ! pressure fraction of cloud top pressure in layer (hPa)
    REAL(jprb)    :: cfraction       ! cloud fraction (0-1) 1 for 100% cloud cover
  END TYPE rttov_profile_aux_s

  !> @internal Auxiliary profile variables
  TYPE rttov_profile_aux
    LOGICAL(jplm) :: on_coef_levels
    TYPE(rttov_profile_aux_s), POINTER :: s(:)

    ! MW CLW absorption
    REAL(jprb), POINTER :: clw(:,:)            ! CLW * layer thickness (g.m^-3.km)

    ! Liebe MW CLW permittivity
    REAL(jprb), POINTER :: debye_prof(:,:,:)   ! Debye terms

    ! Rosenkranz MW CLW permittivity
    REAL(jprb),    POINTER :: ros_eps_s(:,:)   ! Static dielectric constant (epsilon(nu=0))
    REAL(jprb),    POINTER :: ros_dr(:,:)      ! Rosenkranz parameter Delta_R
    REAL(jprb),    POINTER :: ros_gammar(:,:)  ! Rosenkranz parameter gamma_R
    REAL(jprb),    POINTER :: ros_db(:,:)      ! Rosenkranz parameter Delta_B
    REAL(jprb),    POINTER :: ros_nu1(:,:)     ! Rosenkranz parameter nu1
    COMPLEX(jprb), POINTER :: ros_z1(:,:)      ! Rosenkranz parameter z1
    COMPLEX(jprb)          :: ros_z2           ! Rosenkranz parameter z2 (this is a constant)
    REAL(jprb),    POINTER :: ros_log_abs_z1_sq(:,:) ! LOG(ABS(z1)**2)
    COMPLEX(jprb), POINTER :: ros_log_z1star_z2(:,:) ! LOG(CONJG(z1)*z2)
    COMPLEX(jprb), POINTER :: ros_div1(:,:)    ! Precalculated complex division
    COMPLEX(jprb), POINTER :: ros_div2(:,:)    ! Precalculated complex division

    ! TKC MW CLW permittivity
    REAL(jprb),    POINTER :: tkc_eps_s(:,:)   ! Static dielectric constant (epsilon(nu=0))
    REAL(jprb),    POINTER :: tkc_delta_1(:,:) ! TKC parameter delta_1
    REAL(jprb),    POINTER :: tkc_delta_2(:,:) ! TKC parameter delta_2
    REAL(jprb),    POINTER :: tkc_tau_1(:,:)   ! TKC parameter tau_1
    REAL(jprb),    POINTER :: tkc_tau_2(:,:)   ! TKC parameter tau_2

    ! Variables used in water cloud effective diameter parameterisations (visible/IR scattering)
    REAL(jprb), POINTER :: clw_dg(:,:)         ! Generalized effective diameter (microns)

    ! Variables used in ice effective diameter parameterisations (visible/IR scattering)
    REAL(jprb), POINTER :: ice_dg(:,:)         ! Generalized effective diameter (microns)
    REAL(jprb), POINTER :: fac1_ice_dg(:,:)    ! Intermediate variables used to compute the
    REAL(jprb), POINTER :: fac2_ice_dg(:,:)    !   generalized diameter
    REAL(jprb), POINTER :: fac3_ice_dg(:,:)

    ! Variables used in aerosol scattering including rel. hum. calculation
    REAL(jprb),    POINTER :: relhum(:,:)      ! Relative humidity (%)
    REAL(jprb),    POINTER :: tave(:,:)        ! Layer average temperature (K)
    REAL(jprb),    POINTER :: wmixave(:,:)     ! Layer average WV (ppmv)
    REAL(jprb),    POINTER :: xpresave(:,:)    ! Layer average pressure (hPa)
    REAL(jprb),    POINTER :: ppv(:,:)         ! Partial pressure of WV (hPa)
    REAL(jprb),    POINTER :: esw(:,:)         ! Saturated vapour pressure liquid (hPa)
    REAL(jprb),    POINTER :: esi(:,:)         ! Saturated vapour pressure ice (hPa)

    ! Variables used in predictor calculations
    REAL(jprb), POINTER :: t_layer(:,:)        ! avg layer temperature (K)
    REAL(jprb), POINTER :: w_layer(:,:)        ! avg layer humidity (ppmv)
    REAL(jprb), POINTER :: o3_layer(:,:)       ! avg layer ozone (ppmv)
    REAL(jprb), POINTER :: co2_layer(:,:)      ! avg layer co2 (ppmv)
    REAL(jprb), POINTER :: dt(:,:)             ! deviation from ref temperature prof (K)
    REAL(jprb), POINTER :: dto(:,:)            ! deviation from ref ozone temperature prof (K)
    REAL(jprb), POINTER :: tr(:,:), tr_r(:,:), tr_sqrt(:,:) ! ratio t / ref_t
    REAL(jprb), POINTER :: wr(:,:), wr_sqrt(:,:), wr_rsqrt(:,:)
    REAL(jprb), POINTER :: or(:,:), or_sqrt(:,:)
    REAL(jprb), POINTER :: tw(:,:), tw_sqrt(:,:), tw_4rt(:,:)
    REAL(jprb), POINTER :: ww(:,:), ww_r(:,:), wwr(:,:), wwr_r(:,:)
    REAL(jprb), POINTER :: ow(:,:), ow_r(:,:), ow_sqrt(:,:), ow_rsqrt(:,:)
    REAL(jprb), POINTER :: co2r(:,:), co2w(:,:), twr(:,:)
    REAL(jprb), POINTER :: sum(:,:)
  END TYPE rttov_profile_aux

  !> @internal Auxiliary profile variables for RTTOV_SCATT
  TYPE rttov_profile_scatt_aux
    REAL(jprb), POINTER :: cfrac(:)        ! horizontal cloud fraction, one value used for all layers (0-1)
    REAL(jprb), POINTER :: ems_bnd(:)      ! surface emissivity for boundary conditions
    REAL(jprb), POINTER :: ref_bnd(:)      ! surface emissivity for boundary conditions
    REAL(jprb), POINTER :: ems_cld(:)      ! surface emissivity taking into account cloud/rain impact on od
    REAL(jprb), POINTER :: ref_cld(:)      ! surface reflectivity taking into account cloud/rain impact on od
    REAL(jprb), POINTER :: dz(:,:)         ! layer depth   [km]
    REAL(jprb), POINTER :: tbd(:,:)        ! (effective) temperature at layer boundary [K]
    REAL(jprb), POINTER :: tsfc(:)         ! (effective) temperature at surface [K]
    REAL(jprb), POINTER :: tcosmic(:)      ! (effective) temperature of cosmic background [K]
    REAL(jprb), POINTER :: clw(:,:)        ! cloud liquid water (g/m3)
    REAL(jprb), POINTER :: ciw(:,:)        ! cloud ice water (g/m3)
    REAL(jprb), POINTER :: totalice(:,:)   ! total ice (g/m3)
    REAL(jprb), POINTER :: rain(:,:)       ! rain (g/m3)
    REAL(jprb), POINTER :: sp(:,:)         ! solid precipitation (g/m3)
    INTEGER(jpim), POINTER :: mclayer(:)   ! upper level cloud layer
    REAL(jprb), POINTER :: delta(:,:)      ! (= ext*dz/coszen)
    REAL(jprb), POINTER :: tau(:,:)        ! optical depths (= exp(-delta))
    REAL(jprb), POINTER :: ext(:,:)        ! extinction coefficient integreated over hydrometeor types
    REAL(jprb), POINTER :: ssa(:,:)        ! single scattering albedo integreated over hydrometeor types
    REAL(jprb), POINTER :: asm(:,:)        ! asymetry parameter integreated over hydrometeor types [-1,1]
    REAL(jprb), POINTER :: lambda(:,:)     ! eddington approx. variable
                                           ! (= sqrt( 3*ext*ext*(1-ssa)*(1-ssa*asm) )
    REAL(jprb), POINTER :: h (:,:)         ! boundary condition variable (= 1.5_JPRB*ext(1-ssa*asm))
    REAL(jprb), POINTER :: b0(:,:)         ! lower level planck function
    REAL(jprb), POINTER :: b1(:,:)         ! planck function gradient
    REAL(jprb), POINTER :: bn(:,:)         ! upper level planck function
    REAL(jprb), POINTER :: btop(:)         ! planck function entering atmosphere (cosmic radiation)
    REAL(jprb), POINTER :: bsfc(:)         ! planck function at surface  
  END TYPE rttov_profile_scatt_aux

  !> @internal Auxillary variables for RTTOV_SCATT (emissivity retrieval - internal)
  TYPE eddington_sfc_type
    LOGICAL              :: lradiance ! True: up and down are radiances, false: brightness temp.
    REAL(jprb), POINTER  :: tau(:)
    REAL(jprb), POINTER  :: up(:)    ! TOA upwelling radiance or TB
    REAL(jprb), POINTER  :: down(:)  ! SFC downwelling radiance or TB
  END TYPE eddington_sfc_type

  !> @internal Path optical depths as predicted or interpolated (unitless)
  TYPE rttov_opdp_path
    REAL(jprb), POINTER :: atm_level(:,:) ! neg optical depth for thermal radiation (levels to space),
                                          ! size (levels, channels)
    REAL(jprb), POINTER :: sun_level_path2(:,:) ! neg optical depth for solar radiation (levels to space) for
                                                ! combined sun-surface-satellite path, size (levels, channels)
  END TYPE rttov_opdp_path

  !> @internal Transmissions and optical depths (unitless)
  TYPE rttov_path_transmission
    REAL(jprb), POINTER  :: tau_surf_p(:,:)         ! Lambertian transmittance from surface (streams,channels)
    REAL(jprb), POINTER  :: tau_surf_p_r(:,:)       ! reciprocal Lambertian transmittance from surface (streams,channels)
    REAL(jprb), POINTER  :: tau_surf(:,:)           ! transmittance from surface (streams,channels)
    REAL(jprb), POINTER  :: tau_surf_r(:,:)         ! reciprocal transmittance from surface (streams,channels)
    REAL(jprb), POINTER  :: tau_level(:,:,:)        ! transmittance from each standard pressure level
                                                    ! (levels,streams,channels)
    REAL(jprb), POINTER  :: tau_level_r(:,:,:)      ! reciprocal transmittance from each standard pressure level
                                                    ! (levels,streams,channels)
    REAL(jprb), POINTER  :: tau_level_p(:,:,:)      ! Lambertian transmittance from each standard pressure level
                                                    ! (levels,streams,channels)
    REAL(jprb), POINTER  :: tau_level_p_r(:,:,:)    ! reciprocal Lambertian transmittance from each standard pressure level
                                                    ! (levels,streams,channels)
    REAL(jprb), POINTER  :: od_singlelayer(:,:,:)   ! single-layer optical depth
    REAL(jprb), POINTER  :: od_singlelayer_r(:,:,:) ! reciprocal single-layer optical depth
    REAL(jprb), POINTER  :: od_sfrac(:,:)           ! optical depth of partial layer above surface
    REAL(jprb), POINTER  :: od_sfrac_r(:,:)         ! reciprocal optical depth of partial layer above surface
    REAL(jprb), POINTER  :: od_frac_ac(:,:)         ! cloud/aerosol optical depth of partial layer above surface
    REAL(jprb), POINTER  :: tau_surf_ac(:,:)        ! cloud/aerosol transmittance from surface

    REAL(jprb), POINTER  :: fac2(:,:,:)             ! Mask for integration calculation: thermal and solar path1
  END TYPE

  !> @internal Auxiliary transmittance variables
  TYPE rttov_transmission_aux
    REAL(jprb), POINTER  :: fac1(:,:,:)    ! Mask for integration calculation
    REAL(jprb), POINTER  :: surf_fac(:,:)  ! Mask for near-surface layer integration calculation

    TYPE(rttov_path_transmission), POINTER :: thermal_path1   ! Thermal transmittances on surface-sensor path
    TYPE(rttov_path_transmission), POINTER :: solar_path2     ! Solar transmittances on combined sun-surface-sensor path
    TYPE(rttov_path_transmission), POINTER :: solar_path1     ! Solar transmittances on surface-sensor path
  END TYPE rttov_transmission_aux

  !> @internal Auxiliary radiance variables
  TYPE rttov_radiance_aux
    ! Auxiliary calculation arrays for RTE integration
    ! Direct model arrays need to be passed to TL AD and K codes
    ! array size is of (nchannels) or (nlevels, nchannels)
    REAL(jprb), POINTER :: air(:,:)               ! Planck emission from atmospheric layers (mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: surfair(:)             ! Planck emission from near-surface layer (mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: skin(:)                ! Planck emission from surface skin (mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: cosmic(:)              ! Planck emission from CMBR (mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: air_t_eff(:,:)         ! Effective (band-corrected) temperature of atmospheric layers (K)
    REAL(jprb), POINTER :: surf_t_eff(:)          ! Effective (band-corrected) temperature of near-surface layer (K)
    REAL(jprb), POINTER :: skin_t_eff(:)          ! Effective (band-corrected) temperature of surface skin (K)
    REAL(jprb), POINTER :: cosmic_t_eff(:)        ! Effective (band-corrected) temperature of CMBR (K)

    ! Radiances (units: mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: up(:,:,:)               ! sum( B * dT )
    REAL(jprb), POINTER :: down(:,:,:)             ! sum ( B / T**2 dT )
    REAL(jprb), POINTER :: up_solar(:,:,:)         ! sum( B * dT )
    REAL(jprb), POINTER :: down_solar(:,:,:)       ! sum ( B / T**2 dT )
    REAL(jprb), POINTER :: meanrad_up(:,:)
    REAL(jprb), POINTER :: meanrad_down(:,:)
    REAL(jprb), POINTER :: meanrad_up_solar(:,:)
    REAL(jprb), POINTER :: meanrad_down_solar(:,:)
    REAL(jprb), POINTER :: down_ref(:,:,:)
    REAL(jprb), POINTER :: down_ref_solar(:,:,:)
    REAL(jprb), POINTER :: cloudy(:,:)

    ! Quantities used for solar single-scattering
    REAL(jprb), POINTER :: FAC1_2(:,:,:)
    REAL(jprb), POINTER :: FAC3_2(:,:)
    REAL(jprb), POINTER :: FAC4_2(:,:,:)
    REAL(jprb), POINTER :: FAC5_2(:,:,:)
    REAL(jprb), POINTER :: FAC6_2(:,:,:)
    REAL(jprb), POINTER :: FAC7_2(:,:)
    REAL(jprb), POINTER :: FAC4_3(:,:)
    REAL(jprb), POINTER :: FAC5_3(:,:)
  END TYPE rttov_radiance_aux

  !> @internal Raytracing variables
  TYPE rttov_raytracing
    REAL(jprb), POINTER :: ltick(:,:)          ! (levels,profiles) Layer thickness (km)
    REAL(jprb), POINTER :: hgpl(:,:)           ! (levels,profiles) Level geopotential height (km)
    REAL(jprb), POINTER :: dmair(:,:)          ! (levels,profiles) Density of moist air (kg/m3)
    REAL(jprb), POINTER :: refractivity(:,:)   ! (levels,profiles) Refractive index of air
    REAL(jprb), POINTER :: r(:,:)              ! (levels,profiles)
    REAL(jprb), POINTER :: r_r(:,:)            ! (levels,profiles)
    REAL(jprb), POINTER :: z_r(:,:)            ! (layers,profiles)
    REAL(jprb), POINTER :: ratoesun(:,:)       ! (layers,profiles)
    REAL(jprb), POINTER :: ratoesat(:,:)       ! (layers,profiles)
    REAL(jprb), POINTER :: zasun(:,:)          ! (layers,profiles) Sine of local angle of sun-surface path
    REAL(jprb), POINTER :: zasat(:,:)          ! (layers,profiles) Sine of local angle of surface-satellite path
    REAL(jprb), POINTER :: int(:,:)            ! (levels,profiles) Integrated layer values for dmair (hPa/(kg/m3))
    REAL(jprb), POINTER :: ztemp(:,:)          ! (levels,profiles)
    REAL(jprb), POINTER :: ppw(:,:)            ! (levels,profiles) Partial pressure of water vapour (hPa)
    REAL(jprb), POINTER :: dispco2(:,:)        ! (levels,profiles) Correction factor for ref. index due to CO2
    REAL(jprb), POINTER :: pathsat(:,:)        ! (layers,profiles) Secant of local angle of surface-satellite path
    REAL(jprb), POINTER :: pathsat_rsqrt(:,:)  ! (layers,profiles) Reciprocal of SQRT of pathsat
    REAL(jprb), POINTER :: pathsat_sqrt(:,:)   ! (layers,profiles) SQRT of pathsat
    REAL(jprb), POINTER :: pathsun(:,:)        ! (layers,profiles) Secant of local angle of sun-surface path
    REAL(jprb), POINTER :: patheff(:,:)        ! (layers,profiles) Sum of pathsat and pathsun
    REAL(jprb), POINTER :: co2_cm(:)           ! (profiles) CO2 thickness at STP
  END TYPE rttov_raytracing

  !> @internal Sea-surface solar BRDF model variables
  TYPE rttov_sunglint_s
    ! Profile-related quantities
    REAL(jprb) :: windsp
    REAL(jprb) :: wangl
    REAL(jprb) :: dazng
    REAL(jprb) :: zensat
    REAL(jprb) :: zensun

    ! Variables for Yoshimori wave facet model
    REAL(jprb) :: gamma_sq
    REAL(jprb) :: gamma_o
    REAL(jprb) :: gamma_p
    REAL(jprb) :: csi
    REAL(jprb) :: alfa
    REAL(jprb) :: omega
    REAL(jprb) :: gammax
    REAL(jprb) :: q_shad
    REAL(jprb) :: a_shad
    REAL(jprb) :: b_shad
    REAL(jprb) :: lambda_a
    REAL(jprb) :: lambda_b
    REAL(jprb) :: c_shad
    REAL(jprb) :: p_prime
    REAL(jprb) :: g_shad
    REAL(jprb) :: fac1
    REAL(jprb) :: pxy_gammaxy
    REAL(jprb) :: glint

    ! Variables for JONSWAP spectrum
    REAL(jprb) :: x_u
    REAL(jprb) :: alfa1
    REAL(jprb) :: omega_m

    ! Variables for Elfouhaily et al spectrum
    REAL(jprb) :: omega_c
  END TYPE rttov_sunglint_s

  !> @internal Sea-surface solar BRDF model variables
  TYPE rttov_sunglint
    TYPE(rttov_sunglint_s), POINTER :: s(:)   ! (nprofiles)

    ! Variables for JONSWAP spectrum
    REAL(jprb), POINTER :: beta(:,:)          ! (nomega,nprofiles)
    REAL(jprb), POINTER :: psi(:,:)           ! (nomega,nprofiles)

    ! Variables for Elfouhaily et al spectrum
    REAL(jprb), POINTER :: c(:,:)             ! (nk,nprofiles)
    REAL(jprb), POINTER :: lpm(:,:)           ! (nk,nprofiles)
    REAL(jprb), POINTER :: gamma_exp(:,:)     ! (nk,nprofiles)
    REAL(jprb), POINTER :: jp(:,:)            ! (nk,nprofiles)
    REAL(jprb), POINTER :: fpexp(:,:)         ! (nk,nprofiles)
    REAL(jprb), POINTER :: fm(:,:)            ! (nk,nprofiles)
    REAL(jprb), POINTER :: dk(:,:)            ! (nk,nprofiles)
    REAL(jprb), POINTER :: sk2(:,:)           ! (0:nk,nprofiles)
  END TYPE rttov_sunglint

  !> @internal Phase function Legendre coefficients for DOM
  TYPE rttov_phasefn_lcoef
    REAL(jprb), POINTER :: legcoef(:,:,:)   ! Phase function Legendre coefficients (0:dom_nstr,0:1,nlay)
  END TYPE rttov_phasefn_lcoef

  !> @internal optical depths and other quantities for aerosols and clouds
  TYPE rttov_scatt_ir_aercld
    REAL(jprb), POINTER :: opdpsca(:,:)     ! nadir scattering optical depth for all particle types (nlay,nchan)
    REAL(jprb), POINTER :: opdpabs(:,:)     ! nadir absorption optical depth for all particle types (nlay,nchan)
    REAL(jprb), POINTER :: opdpscabpr(:,:)  ! nadir scattering optical depth weighted by bpr for all particle types (nlay,nchan)
    REAL(jprb), POINTER :: opdp(:,:)        ! nadir total layer optical depth (nlay,nchan)
    REAL(jprb), POINTER :: opdpsun(:,:)     ! nadir total layer optical depth on solar effective path (nlay,nchan)
    REAL(jprb), POINTER :: phintup(:,:,:)   ! interpolated phase functions for upward scattering (nparticles,nlay,nchan)
    REAL(jprb), POINTER :: phintdo(:,:,:)   ! interpolated phase functions for downward scattering (nparticles,nlay,nchan)
    REAL(jprb), POINTER :: phtotup(:,:)     ! combined interpolated phase function for upward scattering (nlay,nchan)
    REAL(jprb), POINTER :: phtotdo(:,:)     ! combined interpolated phase function for downward scattering (nlay,nchan)
    REAL(jprb), POINTER :: partsca(:,:,:)   ! scattering coefficient for each particle type (km-1) (nparticles,nlay,nchan)
    REAL(jprb), POINTER :: sca(:,:)         ! total layer scattering coefficient (km-1) (nlay,nchan)
    REAL(jprb), POINTER :: partbpr(:,:,:)   ! bpr value for each particle type (nparticles,nlay,nchan)
  END TYPE rttov_scatt_ir_aercld

  !> @internal Visible/IR scattering optical depths and related data
  TYPE rttov_transmission_scatt_ir
    TYPE(rttov_scatt_ir_aercld), POINTER :: aer
    TYPE(rttov_scatt_ir_aercld), POINTER :: cld
    REAL(jprb), POINTER :: phup(:,:,:)               ! aer/aer+cld interp phasefn for upward scattering (0:1,nlay,nchan)
    REAL(jprb), POINTER :: phdo(:,:,:)               ! aer/aer+cld interp phasefn for downward scattering (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpext(:,:,:)            ! aer/aer+cld extinction optical depth (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpabs(:,:,:)            ! nadir aer/aer+cld absorption optical depth (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpsca(:,:,:)            ! nadir aer/aer+cld scattering optical depth (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpac(:,:,:)             ! aer/aer+cld accumulated optical depth (nlev,0:nstream,nchan)
    REAL(jprb), POINTER :: opdpacl(:,:,:)            ! aer/aer+cld layer optical depth (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpacsun(:,:,:)          ! aer/aer+cld accumulated optical depth (solar path) (nlev,0:nstream,nchan)
    REAL(jprb), POINTER :: opdpaclsun(:,:,:)         ! aer/aer+cld layer optical depth (solar path) (0:1,nlay,nchan)
    REAL(jprb), POINTER :: ssa_solar(:,:,:)          ! aer/aer+cld single-scattering albedo (solar path) (0:1,nlay,nchan)
    REAL(jprb), POINTER :: ssa_thermal(:,:,:)        ! aer/aer+cld single-scattering albedo (0:1,nlay,nchan)
    REAL(jprb), POINTER :: layerod_solar(:,:,:)      ! nadir aer/aer+cld layer optical depths w/ gas (solar path) (0:1,nlay,nchan)
    REAL(jprb), POINTER :: layerod_thermal(:,:,:)    ! nadir aer/aer+cld layer optical depths w/ gas (0:1,nlay,nchan)
    TYPE(rttov_phasefn_lcoef), POINTER :: phasefn(:) ! phase function Legendre coefficients (nchan)
  END TYPE rttov_transmission_scatt_ir

  !> @internal PC coefs for each predictor channel set
  TYPE rttov_coef_pccomp1
    INTEGER(jpim)           :: fmv_pc_npred        ! Number of predictors in the regression set
    INTEGER(jpim), POINTER  :: predictindex(:)     ! Predictors channel indices
    REAL(jprb)   , POINTER  :: coefficients(:,:)   ! Regression coefficients
    REAL(jprb)   , POINTER  :: coefficients_t(:,:) ! Regression coefficients transposed
  END TYPE rttov_coef_pccomp1

  !> @internal PC eigenvectors
  TYPE rttov_coef_pccomp2
    REAL(jprb)   , POINTER  :: eigenvectors  (:,:)   ! Eigenvectors
    REAL(jprb)   , POINTER  :: eigenvectors_t(:,:)   ! Transposed Eigenvectors
  END TYPE rttov_coef_pccomp2

  !> @internal PC-RTTOV coefs
  TYPE rttov_coef_pccomp
    INTEGER(jpim)           :: fmv_pc_comp_pc          ! PC file version number
    INTEGER(jpim)           :: fmv_pc_cld              ! File trained for cloudy simulations
    INTEGER(jpim)           :: fmv_pc_aer              ! File trained for aerosol simulations
    INTEGER(jpim)           :: fmv_pc_naer_types       ! Number of aerosol types for aerosol simulations
    INTEGER(jpim)           :: fmv_pc_nlte             ! File trained for NLTE simulations
    INTEGER(jpim)           :: fmv_pc_msets            ! Maximum number of regression sets
    INTEGER(jpim)           :: fmv_pc_bands            ! Number of bands
    INTEGER(jpim)           :: fmv_pc_mnum             ! Maximum number of eigenvectors
    INTEGER(jpim)           :: fmv_pc_mchn             ! Maximum number of channels
    INTEGER(jpim)           :: fmv_pc_nchn             ! Number of channels
    INTEGER(jpim)           :: fmv_pc_nchn_noise       ! Number of channels for which instrument noise is available
    INTEGER(jpim)           :: fmv_pc_nche             ! Number of channels for which emissisity coefs are available
    INTEGER(jpim)           :: fmv_pc_gas              ! Number of gases for which a reference profile is given
    INTEGER(jpim)           :: fmv_pc_gas_lim          ! Number of gases for which min/max limit profiles are given
    INTEGER(jpim), POINTER  :: fmv_pc_sets   (:)       ! Number of regression sets in each band
    INTEGER(jpim), POINTER  :: emiss_chn     (:)       ! Number of channels for which emissivity coefficients are
    REAL   (jprb), POINTER  :: emiss_c1      (:)       ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c2      (:)       ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c3      (:)       ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c4      (:)       ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c5      (:)       ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c6      (:)       ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c7      (:)       ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c8      (:)       ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c9      (:)       ! Emissivity coefficient
    INTEGER(jpim)           :: fmv_pc_nlev             ! Number of reference profile levels
    REAL(jprb), POINTER     :: ref_pc_prfl_p (:)       ! pressure  (hPa)       (levels)
    REAL(jprb), POINTER     :: ref_pc_prfl_mr(:,:)     ! mixing ratio (ppmv)   (levels)
    REAL(jprb), POINTER     :: lim_pc_prfl_tmin(:)     ! Profile limit :temperature (K)
    REAL(jprb), POINTER     :: lim_pc_prfl_tmax(:)     ! Profile limit :temperature (K)
    REAL(jprb), POINTER     :: lim_pc_prfl_qmin(:)     ! Profile limit :water vapour (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_qmax(:)     ! Profile limit :water vapour (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_ozmin(:)    ! Profile limit :ozone (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_ozmax(:)    ! Profile limit :ozone (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_gasmin(:,:) ! Profile limit :additional gases (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_gasmax(:,:) ! Profile limit :additional gases (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_aermin(:,:) ! Profile limit :aerosols (cm^-3) (layers)
    REAL(jprb), POINTER     :: lim_pc_prfl_aermax(:,:) ! Profile limit :aerosols (cm^-3) (layers)
    REAL(jprb)              :: lim_pc_prfl_pmin        ! Surface pressure (hPa)
    REAL(jprb)              :: lim_pc_prfl_pmax        ! Surface pressure (hPa)
    REAL(jprb)              :: lim_pc_prfl_tsmin       ! Surface temperature (K)
    REAL(jprb)              :: lim_pc_prfl_tsmax       ! Surface temperature (K)
    REAL(jprb)              :: lim_pc_prfl_skmin       ! Skin temperature (K)
    REAL(jprb)              :: lim_pc_prfl_skmax       ! Skin temperature (K)
    REAL(jprb)              :: lim_pc_prfl_wsmin       ! 10m wind speed (m/s)
    REAL(jprb)              :: lim_pc_prfl_wsmax       ! 10m wind speed (m/s)
    REAL(jprb), POINTER     :: co2_pc_ref(:)           ! Fixed co2 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: n2o_pc_ref(:)           ! Fixed n2o profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co_pc_ref(:)            ! Fixed co  profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: ch4_pc_ref(:)           ! Fixed ch4 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co2_pc_min(:)           ! Fixed co2 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: n2o_pc_min(:)           ! Fixed n2o profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co_pc_min(:)            ! Fixed co  profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: ch4_pc_min(:)           ! Fixed ch4 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co2_pc_max(:)           ! Fixed co2 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: n2o_pc_max(:)           ! Fixed n2o profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co_pc_max(:)            ! Fixed co  profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: ch4_pc_max(:)           ! Fixed ch4 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: noise_in(:)             ! Noise values for the channels whose radiances are
                                                       ! reconstructed using principal components
    REAL(jprb), POINTER     :: noise(:)                ! Noise values for the channels whose radiances are
                                                       ! used as predictors in the computation of principal components
    REAL(jprb), POINTER     :: noise_r(:)              ! Reciprocal noise
    INTEGER(jpim), POINTER  :: ff_ori_chn_in(:)
    REAL(jprb),    POINTER  :: ff_cwn_in(:)            ! central wave number of reconstructed radiances (cm-1)
    REAL(jprb),    POINTER  :: ff_bco_in (:)           ! band correction offset (K)
    REAL(jprb),    POINTER  :: ff_bcs_in (:)           ! band correction slope (K/K)
    REAL(jprb),    POINTER  :: planck1_in(:)           ! C1 * Nu**3 (mW/(m2*sr*cm-1))
    REAL(jprb),    POINTER  :: planck2_in(:)           ! C2 * Nu (K)
    TYPE(rttov_coef_pccomp1), POINTER:: pcreg(:,:)
    TYPE(rttov_coef_pccomp2), POINTER:: eigen(:)
  END TYPE rttov_coef_pccomp

  !> @internal IR scattering coefs (aer=aerosol, wcl=water cloud, icl=ice cloud)
  TYPE rttov_coef_scatt_ir
    INTEGER(jpim)             :: fmv_aer_chn          ! Number of channels for which optical parameters are stored
    INTEGER(jpim)             :: fmv_wcl_chn
    INTEGER(jpim)             :: fmv_wcldeff_chn      ! Cloud liquid water Deff scheme
    INTEGER(jpim)             :: fmv_icl_chn
    INTEGER(jpim)             :: fmv_aer_pha_chn      ! Number of channels for which phase functions are stored
    INTEGER(jpim)             :: fmv_wcl_pha_chn
    INTEGER(jpim)             :: fmv_wcldeff_pha_chn
    INTEGER(jpim)             :: fmv_icl_pha_chn
    INTEGER(jpim)             :: fmv_aer_comp         ! Number of particle types
    INTEGER(jpim)             :: fmv_wcl_comp
    INTEGER(jpim)             :: fmv_aer_maxnmom      ! Maximum number of Legendre moments over all phase functions
    INTEGER(jpim)             :: fmv_wcl_maxnmom
    INTEGER(jpim)             :: fmv_wcldeff_maxnmom
    INTEGER(jpim)             :: fmv_icl_maxnmom
    INTEGER(jpim)             :: aer_nphangle         ! Size of angle grid for phase functions
    INTEGER(jpim)             :: wcl_nphangle
    INTEGER(jpim)             :: wcldeff_nphangle
    INTEGER(jpim)             :: icl_nphangle
    INTEGER(jpim)             :: fmv_wcldeff_ndeff    ! Size of water cloud effective diameter grid
    INTEGER(jpim)             :: fmv_icl_ndeff        ! Size of ice cloud effective diameter grid
    CHARACTER(LEN=4), POINTER :: fmv_aer_comp_name(:) ! Particle type names
    CHARACTER(LEN=4), POINTER :: fmv_wcl_comp_name(:)
    INTEGER(jpim), POINTER    :: fmv_aer_rh(:)        ! Number of relative humidity values for each particle type
    INTEGER(jpim), POINTER    :: fmv_wcl_rh(:)
    REAL(jprb),    POINTER    :: fmv_aer_rh_val(:)    ! Relative humidity values for each particle type (%)
    REAL(jprb),    POINTER    :: fmv_wcl_rh_val(:)
    REAL(jprb),    POINTER    :: fmv_wcldeff_deff(:)  ! Water cloud effective diameter values (microns)
    REAL(jprb),    POINTER    :: fmv_icl_deff(:)      ! Ice cloud effective diameter values (microns)
    INTEGER(jpim), POINTER    :: aer_pha_chanlist(:)  ! List of channels for which phase functions are stored
    INTEGER(jpim), POINTER    :: wcl_pha_chanlist(:)
    INTEGER(jpim), POINTER    :: wcldeff_pha_chanlist(:)
    INTEGER(jpim), POINTER    :: icl_pha_chanlist(:)
    INTEGER(jpim), POINTER    :: aer_pha_index(:)     ! Channel indices for stored phase fns
    INTEGER(jpim), POINTER    :: wcl_pha_index(:)
    INTEGER(jpim), POINTER    :: wcldeff_pha_index(:)
    INTEGER(jpim), POINTER    :: icl_pha_index(:)
    REAL(jprb),    POINTER    :: aer_phangle(:)       ! Angles over which phase fns are defined
    REAL(jprb),    POINTER    :: wcl_phangle(:)
    REAL(jprb),    POINTER    :: wcldeff_phangle(:)
    REAL(jprb),    POINTER    :: icl_phangle(:)
    TYPE(rttov_phasefn_int)   :: aer_phfn_int         ! Auxiliary data for rapidly interpolating phase functions
    TYPE(rttov_phasefn_int)   :: wcl_phfn_int
    TYPE(rttov_phasefn_int)   :: wcldeff_phfn_int
    TYPE(rttov_phasefn_int)   :: icl_phfn_int
    REAL(jprb),    POINTER    :: abs(:,:)             ! Absorption coefficients normalised to 1 part.cm-3 (km-1/part.cm-3)
    REAL(jprb),    POINTER    :: sca(:,:)             ! Scattering coefficients normalised to 1 part.cm-3 (km-1/part.cm-3)
    REAL(jprb),    POINTER    :: bpr(:,:)             ! Backscattering parameters for Chou-scaling
    INTEGER(jpim), POINTER    :: nmom(:,:)            ! Number of Legendre moments for each phase function
    REAL(jprb),    POINTER    :: legcoef(:,:,:)       ! Legendre coefficients for each phase function
    REAL(jprb),    POINTER    :: pha(:,:,:)           ! Phase functions
    REAL(jprb),    POINTER    :: confac(:)            ! Cloud liquid water unit conversion factors (part.cm-3/g.m-3)
    REAL(jprb),    POINTER    :: aer_mmr2nd(:)        ! Aerosol unit conversion factors (g.m-3/part.cm-3)
  END TYPE rttov_coef_scatt_ir

  !> @internal Baran scheme variables
  TYPE rttov_coef_optpcld
    ! interpolation factors for frequency
    INTEGER(jpim), POINTER  :: iwn(:)
    INTEGER(jpim), POINTER  :: jwn(:)
    REAL(jprb),    POINTER  :: dx_dwn(:)

    REAL(jprb),    POINTER  :: q(:), w(:)     ! Gaussian quadrature to calculate phase fn Leg. coefs
    TYPE(rttov_phasefn_int) :: phfn_int       ! Phase fn interpolation data
  END TYPE rttov_coef_optpcld

  !> @internal IR scattering optical properties
  TYPE rttov_optpar_ir
    TYPE(rttov_coef_scatt_ir), POINTER :: optpaer(:)
    TYPE(rttov_coef_scatt_ir), POINTER :: optpwcl(:)
    TYPE(rttov_coef_scatt_ir), POINTER :: optpwcldeff
    TYPE(rttov_coef_scatt_ir), POINTER :: optpicl
    TYPE(rttov_coef_optpcld),  POINTER :: optpiclbaran
  END TYPE rttov_optpar_ir

  !> @internal MFASIS LUT data for one channel
  TYPE rttov_mfasis_lut
    INTEGER(jpim)       :: nluts       ! number of LUTs for this channel
    REAL(jprb), POINTER :: q(:)        ! water vapour value for each LUT (nluts)
    REAL(jprb), POINTER :: data(:,:)   ! LUT(s) for a single channel (all dimensions concatenated per LUT) (:,nluts)
  END TYPE rttov_mfasis_lut

  !> @internal MFASIS LUT axis definition
  TYPE rttov_mfasis_axis
    CHARACTER(LEN=32)   :: name        ! axis variable name
    INTEGER(jpim)       :: dim_type    ! type of dimension (interp. method):
                                       !   0 : theta+ Fourier index (k)
                                       !   1 : theta- Fourier index (l)
                                       !   2 : Albedo (no interp., use Jonkheid 2012)
                                       !   3 : tau-like (interpolate in logarithm)
                                       !   4 : Reff-like (interpolate linearly)
                                       !   5 : scattering angle (constr. lin. interp.)
    INTEGER(jpim)       :: nvalues
    REAL(jprb), POINTER :: values(:)   ! axis values (nvalues)
  END TYPE rttov_mfasis_axis

  !> @internal MFASIS cloud or aerosol LUT structure
  TYPE rttov_coef_mfasis
    INTEGER(jpim) :: file_type                      ! 1=>clouds, 2=>aerosols
    INTEGER(jpim) :: version                        ! version number, may be useful for handling
                                                    !   updates to LUT format later and checking
                                                    !   compatibility with code
    CHARACTER(LEN=132) :: readme_lut(100)           ! description of LUT creation, etc

    INTEGER(jpim)          :: ndims                 ! number of LUT dimensions
    TYPE(rttov_mfasis_axis), POINTER :: lut_axes(:) ! table axis values (ndims)

    INTEGER(jpim)          :: nparticles            ! number of particle types in the table (for
                                                    !   clouds this is 2, water and ice cloud; for
                                                    !   aerosol this has yet to be decided)
    INTEGER(jpim), POINTER :: aer_types(:)          ! the indices of the RTTOV aerosol types used
                                                    !   in training aerosol LUTs (nparticles)

    INTEGER(jpim)          :: clw_scheme            ! liquid cloud scheme used for training LUT
    INTEGER(jpim)          :: ice_scheme            ! ice cloud scheme used for training LUT

    INTEGER(jpim)          :: nchannels             ! number of supported channels
    INTEGER(jpim), POINTER :: channel_list(:)       ! supported channel numbers (nchannels)
    INTEGER(jpim)          :: nchannels_coef        ! number of channels in the associated rtcoef file
    INTEGER(jpim), POINTER :: channel_lut_index(:)  ! index of each supported channel in lut(:)  (nchannels_coef)
    TYPE(rttov_mfasis_lut), POINTER :: lut(:)       ! LUT(s) for each supported channel (nchannels)
  END TYPE rttov_coef_mfasis

  !> @internal HT-FRTC scheme coefficients (Optical properties + PC related)
  TYPE rttov_coef_htfrtc
    INTEGER(jpim)          :: opt_prop_type
    INTEGER(jpim)          :: n_f
    REAL(jprb), POINTER    :: freq(:)
    INTEGER(jpim)          :: n_gas_l
    INTEGER(jpim), POINTER :: gasid_l(:)
    INTEGER(jpim)          :: n_gas_nl
    INTEGER(jpim), POINTER :: gasid_nl(:)
    INTEGER(jpim)          :: n_p
    REAL(jprb), POINTER    :: p(:)
    INTEGER(jpim)          :: n_val_l
    INTEGER(jpim)          :: n_val_nl
    REAL(jprb), POINTER    :: val_l(:)
    REAL(jprb), POINTER    :: val_nl(:)
    INTEGER(jpim)          :: n_t
    REAL(jprb), POINTER    :: val_t(:)
    REAL(jprb), POINTER    :: val_mf_l(:)
    INTEGER(jpim)          :: n_mf_nl
    REAL(jprb), POINTER    :: val_mf_nl(:,:)
    INTEGER(jpim)          :: n_b
    REAL(jprb), POINTER    :: val_b(:)
    INTEGER(jpim)          :: n_lt
    REAL(jprb), POINTER    :: val_lt(:)
    REAL(jprb), POINTER    :: coef_l(:,:,:,:)
    REAL(jprb), POINTER    :: coef_nl(:,:,:,:)
    REAL(jprb), POINTER    :: coef_b(:,:)
    REAL(jprb), POINTER    :: coef_lt(:)
    INTEGER(jpim)          :: n_ssemp
    REAL(jprb), POINTER    :: coef_ssemp(:,:)
    INTEGER(jpim)          :: n_pc
    REAL(jpim), POINTER    :: coef_pdt(:,:)
    REAL(jpim), POINTER    :: val_mean(:)
    REAL(jpim), POINTER    :: val_norm(:)
    INTEGER(jpim)          :: n_ch
    REAL(jprb), POINTER    :: sensor_freq(:)
    REAL(jprb), POINTER    :: ch_mean(:)
    REAL(jpim), POINTER    :: pc(:,:)
    REAL(jprb), POINTER    :: coef_aux(:)
  END TYPE rttov_coef_htfrtc

  !> @internal RTTOV coefs
  TYPE rttov_coefs
    LOGICAL(jplm)              :: initialised = .FALSE.
    TYPE(rttov_coef)           :: coef
    TYPE(rttov_coef_scatt_ir)  :: coef_scatt_ir
    TYPE(rttov_optpar_ir)      :: optp
    TYPE(rttov_coef_pccomp)    :: coef_pccomp
    TYPE(rttov_coef_mfasis)    :: coef_mfasis_cld
    TYPE(rttov_coef_mfasis)    :: coef_mfasis_aer
    TYPE(rttov_coef_htfrtc)    :: coef_htfrtc
  END TYPE rttov_coefs

  !> @internal IR cloud stream variables
  TYPE rttov_ircld
    INTEGER(jpim), POINTER  :: nstream(:)
    INTEGER(jpim), POINTER  :: nstreamref(:)
    INTEGER(jpim), POINTER  :: iloop(:)
    INTEGER(jpim), POINTER  :: icount(:)
    INTEGER(jpim), POINTER  :: icounstr(:)
    INTEGER(jpim), POINTER  :: icount1(:)
    REAL(jprb)   , POINTER  :: xstrclr(:)
    INTEGER(jpim), POINTER  :: icldarr   (:,:,:)
    REAL(jprb)   , POINTER  :: xstrref1  (:,:,:)
    REAL(jprb)   , POINTER  :: xstrref2  (:,:,:)
    INTEGER(jpim), POINTER  :: indexstr  (:,:)
    INTEGER(jpim), POINTER  :: icount1ref(:,:)
    INTEGER(jpim), POINTER  :: iloopin   (:,:)
    INTEGER(jpim), POINTER  :: iflag     (:,:)
    REAL(jprb)   , POINTER  :: xstr      (:,:)
    REAL(jprb)   , POINTER  :: xstrminref(:,:)
    REAL(jprb)   , POINTER  :: xstrref   (:,:)
    REAL(jprb)   , POINTER  :: cldcfr    (:,:)
    REAL(jprb)   , POINTER  :: maxcov    (:,:)
    REAL(jprb)   , POINTER  :: xstrmax   (:,:)
    REAL(jprb)   , POINTER  :: xstrmin   (:,:)
    REAL(jprb)   , POINTER  :: a         (:,:)
    REAL(jprb)   , POINTER  :: ntotref   (:,:)
    LOGICAL(jplm), POINTER  :: flag      (:,:)
  END TYPE rttov_ircld

  !> @internal Profiles input to Discrete Ordinates algorithm
  TYPE rttov_profile_dom
    INTEGER(jpim)          :: nlayers       ! Number of layers in DOM profile
    LOGICAL(jplm)          :: surface       ! Flag to indicate if surface is visible
    REAL(jprb),    POINTER :: layerod(:)    ! Layer total optical depth (nlayers)
    INTEGER(jpim), POINTER :: laymap(:)     ! Mapping from DOM layer to user layer
  END TYPE

  !> @internal Direct model internal state of Discrete Ordinates algorithm
  TYPE rttov_dom_state
    ! Anything of size 0:naz is 0:0 for thermal channels, 0:nstr-1 for solar

    ! Solar only
    INTEGER(jpim), POINTER :: nazloops(:)      ! 0:nstreams
    REAL(jprb),    POINTER :: z(:,:,:,:)       ! nstr,0:profiles(1)%nlayers,0:1,0:nstr-1

    ! Thermal only
    REAL(jprb),    POINTER :: y0(:,:,:)        ! nstr,nlayers,0:nstreams
    REAL(jprb),    POINTER :: y1(:,:,:)        ! nstr,nlayers,0:nstreams

    ! Thermal and solar
    REAL(jprb),    POINTER :: thisrad(:)       ! 0:nstreams
    REAL(jprb),    POINTER :: radsurfup(:)     ! 0:nstreams
    REAL(jprb),    POINTER :: radsurfup_sum(:) ! 0:nstreams
    REAL(jprb),    POINTER :: x(:,:,:)         ! nstr*maxnlayers, 0:naz, 0:nstreams
    REAL(jprb),    POINTER :: xp(:,:,:,:,:)    ! nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1
    REAL(jprb),    POINTER :: eval(:,:,:,:)    ! nstr/2,profiles(1)%nlayers,0:1,0:naz
  END TYPE

  !> @internal Direct model internal state of MFASIS
  TYPE rttov_mfasis_refl
    REAL(jprb)             :: refl             ! Reflectance computed by MFASIS
    REAL(jprb),    POINTER :: refl_lin_coef(:) ! Sensitivity of reflectance (matrix elements of linear computations)
  END TYPE rttov_mfasis_refl

  !> @internal RTTOV internal state
  TYPE rttov_traj
!
! Hold RTTOV trajectory; these variables have counterparts in TL, AD, K,
! and/or their dimensions are known before running RTTOV (nlevels, nprofiles, nchannels)
! it is possible to allocate these variables from outside RTTOV
!
    TYPE(rttov_profile),      POINTER :: profiles_coef(:)
    TYPE(rttov_profile),      POINTER :: profiles_int(:)
    TYPE(rttov_predictors)            :: predictors
    TYPE(rttov_profile_aux)           :: aux_prof
    TYPE(rttov_profile_aux)           :: aux_prof_coef
    TYPE(rttov_raytracing)            :: raytracing
    TYPE(rttov_raytracing)            :: raytracing_coef
    TYPE(rttov_opdp_path)             :: opdp_path
    TYPE(rttov_opdp_path)             :: opdp_path_coef

    REAL(jprb),               POINTER :: diffuse_refl(:) ! Surface refl for downwelling radiation (nchanprof)
    REAL(jprb),               POINTER :: fresnrefl(:)    ! Fresnel reflection coefficients (nchanprof)
    TYPE(rttov_sunglint)              :: sunglint
    TYPE(rttov_ircld)                 :: ircld
    TYPE(rttov_transmission_scatt_ir) :: transmission_scatt_ir

    TYPE(rttov_radiance_aux)          :: auxrad

    TYPE(rttov_coefs),        POINTER :: coefs
    INTEGER(jpim)                     :: nchanprof
    INTEGER(jpim)                     :: nlevels
    INTEGER(jpim)                     :: nlayers
    TYPE(rttov_options)               :: opts
  END TYPE rttov_traj

  !> @internal RTTOV internal state
  TYPE rttov_traj_dyn
!
! Hold RTTOV trajectory; these variables have counterparts in TL, AD, K,
! but their dimensions are known when RTTOV starts running (nstreams)
!
    LOGICAL(jplm)                     :: from_tladk = .FALSE.  ! True if this direct model traj_dyn came from the TL/AD/K
    INTEGER(jpim)                     :: nstreams = -1  ! This initialisation used to determine alloc status
    TYPE(rttov_radiance_aux)          :: auxrad_stream
    TYPE(rttov_transmission_scatt_ir) :: transmission_scatt_ir_dyn
    TYPE(rttov_transmission_aux)      :: transmission_aux
    TYPE(rttov_profile_dom), POINTER  :: profiles_dom_thermal(:,:)
    TYPE(rttov_profile_dom), POINTER  :: profiles_dom_solar(:,:)
    TYPE(rttov_dom_state),   POINTER  :: dom_state_thermal(:)
    TYPE(rttov_dom_state),   POINTER  :: dom_state_solar(:)
    TYPE(rttov_mfasis_refl), POINTER  :: mfasis_refl(:,:)
  END TYPE rttov_traj_dyn

  !> @internal RTTOV internal state (optical depths, transmittances)
  TYPE rttov_path_traj_trans
    ! Structure to hold optical depth and transmittance data
    ! within static trajectory.
    REAL(jprb),     POINTER :: tau_ref       (:,:)
    REAL(jprb),     POINTER :: tau_ref_surf  (:)
    REAL(jprb),     POINTER :: tau_level     (:,:)! sat to level transmittance
    REAL(jprb),     POINTER :: tau_surf      (:)
    REAL(jprb),     POINTER :: od_level      (:,:)! sat to level optical depth
    REAL(jprb),     POINTER :: opdp_ref_coef (:,:)! layer optical depth before threshold

    REAL(jprb),     POINTER :: od_singlelayer(:,:)! single layer optical depth
    REAL(jprb),     POINTER :: od_frac       (:)
  END TYPE rttov_path_traj_trans

  !> @internal RTTOV internal state (direct model only)
  TYPE rttov_traj_sta
!
! Hold RTTOV trajectory; these variables do not have counterparts in TL, AD, K
!
    LOGICAL(jplm)                   :: do_opdep_calc        ! flag to indicate RTTOV gas optical depth calculation required
    LOGICAL(jplm)                   :: do_mfasis            ! flag to indicate MFASIS is being used

    LOGICAL(jplm),  POINTER         :: thermal(:)           ! switch for thermal calculations (nchanprof)
    LOGICAL(jplm),  POINTER         :: solar(:)             ! switch for solar calculations (nchanprof)
    LOGICAL(jplm)                   :: dothermal            ! flag to indicate thermal calculations required
    LOGICAL(jplm)                   :: dosolar              ! flag to indicate solar calculations required

    INTEGER(jpim)                   :: dom_nstreams         ! number of streams for DOM
    LOGICAL(jplm)                   :: plane_parallel       ! switch for plane-parallel geometry (no atmospheric curvature)

    REAL(jprb), POINTER             :: solar_spec_esd(:)    ! Solar spectrum adjusted for esd (nchanprof) (mW/m2/cm-1)
    REAL(jprb), POINTER             :: refl_norm(:)         ! Normalisation factor for solar surface reflectance

    LOGICAL(jplm), POINTER          :: do_lambertian(:)
    TYPE(rttov_path_traj_trans), POINTER :: thermal_path1
    TYPE(rttov_path_traj_trans), POINTER :: solar_path2
    TYPE(rttov_path_traj_trans), POINTER :: solar_path1
    TYPE(rttov_geometry), POINTER        :: angles(:)
    TYPE(rttov_geometry), POINTER        :: angles_coef(:)
    TYPE(rttov_profile),  POINTER        :: profiles_coef_ref(:)
    TYPE(rttov_chanprof), POINTER        :: chanprof_in(:)
    TYPE(rttov_chanprof), POINTER        :: chanprof_pc(:)
    REAL(jprb), POINTER                  :: pc_aer_ref(:,:,:) ! For PC-RTTOV aerosol simulations
    REAL(jprb), POINTER                  :: pc_aer_min(:,:,:) !
    REAL(jprb), POINTER                  :: pc_aer_max(:,:,:) !
  END TYPE rttov_traj_sta

  !> @internal Used for coef testing
  TYPE rttov_lbl_check
    REAL(jprb), POINTER :: atm_layer(:,:)
    REAL(jprb), POINTER :: atm_layer_path2(:,:)
    LOGICAL(jplm)       :: plane_geometry
  END TYPE rttov_lbl_check

END MODULE rttov_types

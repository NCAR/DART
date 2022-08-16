.. _obs_def_rttov_mod:

MODULE ``obs_def_rttov_mod``
============================

Overview
--------

DART RTTOV observation module, including the observation operators for the two primary 
RTTOV-observation types -- visible/infrared radiances and microwave 
radiances/brightness temperatures.

The obs_def_rttov_mod.f90 module acts as a pass-through for RTTOV version 12.3. For more information,
see `the RTTOV site <https://www.nwpsaf.eu/site/software/rttov/documentation/>`__.

For RTTOV v13 use the obs_def_rttov13_mod.f90 module contributed by Lukas Kugler
of the University of Vienna.

DART supports both RTTOV-direct for visible/infrared/microwave as well as RTTOV-scatt 
for microwave computations. The code, in principle, supports all features of version 12.3 
as a pass-through from the model to RTTOV, includes aerosols, trace gases, clouds, and 
atmospheric variables. The code also includes directly specifying scattering properties.

However, a model may not have all of the variables necessary for these functions 
depending on your model's setup.  For example, DART can use any of the RTTOV clw or ice 
schemes, but the WRF model is not directly compatible with the IR default cloud 
classification of marine/continental stratus/cumulus clean/dirty. We also offer a simple
classification based on maximum vertical velocity in the column and land type, but due to 
lack of aerosol information, WRF/DART cannot differentiate between clean and dirty cumulus. 
This may have some impact on the forward calculations - but in experience the difference 
in cloud phase (ice versus water) makes a much larger difference.  Trace gases and aerosols 
may be important for actual observation system experiments using visible/infrared; this may
depend on the precise frequencies you wish to use.

Although a model may not have the necessary inputs by itself,
the defaults in RTTOV based on climatology can be used.
The impact on the quality of the results should be investigated.

Known issues:

-  DART does not yet provide any type of bias correction
-  Cross-channel error correlations are not yet supported. A principal component approach has been discussed. For now,
   the best bet is to use a subset of channels that are nearly independent of one another.
-  Vertical localization will need to be tuned. Turning off vertical localization may work well if you have a large
   number of ensemble members. Using the maximum peak of the weighting function or the cloud-top may be appropriate.
   There are also other potential approaches being investigated.

| Author and Contact information:

-  DART Code: Jeff Steward
-  Original DART/RTTOV work: Nancy Collins, Johnny Hendricks

Backward compatibility note
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Other modules used
------------------

::

   types_mod
   utilities_mod
   location_mod (threed_sphere)
   assim_model_mod
   obs_def_utilitie_mod
   ensemble_manager_mod
   utilities_mod
   parkind1 (from RTTOV)
   rttov_types (from RTTOV)
   obs_kind_mod

Public interfaces
-----------------

=============================== ========================
*use obs_def_rttov_mod, only :* set_visir_metadata
\                               set_mw_metadata
\                               get_expected_radiance
\                               get_rttov_option_logical
=============================== ========================

Namelist interface ``&obs_def_rttov_mod_nml`` is read from file ``input.nml``.

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, & platform_id, sat_id, sensor_id, channel,
   specularity)*
   ::

      integer,  intent(out) :: key
      real(r8), intent(in)  :: sat_az
      real(r8), intent(in)  :: sat_ze
      real(r8), intent(in)  :: sun_az
      real(r8), intent(in)  :: sun_ze
      integer,  intent(in)  :: platform_id, sat_id, sensor_id, channel
      real(r8), intent(in)  :: specularity

.. container:: indent1

   Visible / infrared observations have several auxillary metadata variables. Other than the key, which is standard DART
   fare, the RTTOV satellite azimuth and satellite zenith angle must be specified. See the RTTOV user guide for more
   information (in particular, see figure 4). If the ``addsolar`` namelist value is set to true, then the solar azimuth
   and solar zenith angles must be specified - again see the RTTOV user guide. In addition to the platform/satellite/
   sensor ID numbers, which are the RTTOV unique identifiers, the channel specifies the chanenl number in the RTTOV
   coefficient file. Finally, if ``do_lambertian`` is true, specularity must be specified here. Again, see the RTTOV
   user guide for more information.

   =============== ================================================================
   ``key``         The DART observation key.
   ``sat_az``      The satellite azimuth angle.
   ``sat_ze``      The satellite zenith angle.
   ``sun_az``      The solar azimuth angle. Only relevant if addsolar is true.
   ``sun_ze``      The solar zenith angle. Only relevant if addsolar is true.
   ``platform_id`` The RTTOV platform ID.
   ``sat_id``      The RTTOV satellite ID.
   ``sensor_id``   The RTTOV sensor ID.
   ``channel``     The RTTOV channel number.
   ``specularity`` The surface specularity. Only relevant if do_lambertian is true.
   =============== ================================================================

| 

.. container:: routine

   *call set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, channel, mag_field, cosbk, fastem_p1,
   fastem_p2, fastem_p3, fastem_p4, fastem_p5)*
   ::

      integer,  intent(out) :: key
      real(r8), intent(in)  :: sat_az
      real(r8), intent(in)  :: sat_ze
      integer,  intent(in)  :: platform_id, sat_id, sensor_id, channel
      real(r8), intent(in)  :: mag_field
      real(r8), intent(in)  :: cosbk
      real(r8), intent(in)  :: fastem_p[1-5]

.. container:: indent1

   Microwave observations have several auxillary metadata variables. Other than the key, which is standard DART fare,
   the RTTOV satellite azimuth and satellite zenith angle must be specified. See the RTTOV user guide for more
   information (in particular, see figure 4). In addition to the platform/satellite/ sensor ID numbers, which are the
   RTTOV unique identifiers, the channel specifies the chanenl number in the RTTOV coefficient file. In addition, if
   ``use_zeeman`` is true, the magnetic field and cosine of the angle between the magnetic field and angle of
   propagation must be specified. See the RTTOV user guide for more information. Finally, the fastem parameters for land
   must be specified here. This may be difficult for observations to set, so default values (see table 21 in the RTTOV
   user guide) can be used until a better solution is devised.

   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``key``           | The DART observation key.                                                                       |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``sat_az``        | The satellite azimuth angle.                                                                    |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``sat_ze``        | The satellite zenith angle.                                                                     |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``platform_id``   | The RTTOV platform ID.                                                                          |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``sat_id``        | The RTTOV satellite ID.                                                                         |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``sensor_id``     | The RTTOV sensor ID.                                                                            |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``channel``       | The RTTOV channel number.                                                                       |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``mag_field``     | The strength of the magnetic field. Only relevant if add_zeeman is true.                        |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``cosbk``         | The cosine of the angle between the magnetic field and direction of EM propagation. Only        |
   |                   | relevant if add_zeeman is true.                                                                 |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | ``fastem_p[1-5]`` | The five parameters used for fastem land/sea ice emissivities. For ocean emissivities, an       |
   |                   | internal model is used based on the value of fastem_version.                                    |
   +-------------------+-------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call get_expected_radiance(obs_kind_ind, state_handle, ens_size, location, key, val, istatus)*
   ::

      integer,             intent(in)  :: obs_kind_ind
      type(ensemble_type), intent(in)  :: state_handle
      integer,             intent(in)  :: ens_size
      type(location_type), intent(in)  :: location
      integer,             intent(in)  :: key
      real(r8),            intent(out) :: val(ens_size)
      integer,             intent(out) :: istatus(ens_size)

.. container:: indent1

   Given a location and the state vector from one of the ensemble members, compute the model-predicted satellite
   observation. This can be either in units of radiance (mW/cm-1/sr/sq.m) or a brightness temperature (in K), depending
   on if this is a visible/infrared observation or a microwave observation.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``obs_kind_ind`` | The index of the observation kind; since many observation kinds are handled by this module, this |
   |                  | can be used to determine precisely which observation kind is being used.                         |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``state_handle`` | The ensemble of model states to be used for the observation operator calculations.               |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``location``     | Location of this observation                                                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``key``          | Unique identifier associated with this satellite observation                                     |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``val``          | The returned observation in units of either radiance or brightness temperature.                  |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``istatus``      | Returned integer status code describing problems with applying forward operator. 0 is a good     |
   |                  | value; any positive value indicates an error; negative values are reserved for internal DART use |
   |                  | only.                                                                                            |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *p = get_rttov_option_logical(field_name)*
   ::

      character(len=*),           intent(in)  :: field_name
      logical,                    result      :: p

.. container:: indent1

   Return the logical value of the RTTOV parameter associated with the field_name.

   ============== =======================================================
   ``field_name`` The name of the RTTOV parameter from the namelist.
   ``p``          The logical return value associated with the parameter.
   ============== =======================================================

| 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_def_rttov_nml
      rttov_sensor_db_file   = 'rttov_sensor_db.csv'
      first_lvl_is_sfc       = .true. 
      mw_clear_sky_only      = .false.
      interp_mode            = 1 
      do_checkinput          = .true.
      apply_reg_limits       = .true.
      verbose                = .true.
      fix_hgpl               = .false.
      do_lambertian          = .false.
      lambertian_fixed_angle = .true.
      rad_down_lin_tau       = .true.
      use_q2m                = .true.
      use_uv10m              = .true.
      use_wfetch             = .false.
      use_water_type         = .false.
      addrefrac              = .false.
      plane_parallel         = .false.
      use_salinity           = .false.
      apply_band_correction  = .true.
      cfrac_data             = .true.
      clw_data               = .true.
      rain_data              = .true.
      ciw_data               = .true.
      snow_data              = .true.
      graupel_data           = .true.
      hail_data              = .false.
      w_data                 = .true.
      clw_scheme             = 1
      clw_cloud_top          = 322.
      fastem_version         = 6
      supply_foam_fraction   = .false.
      use_totalice           = .true.
      use_zeeman             = .false.
      cc_threshold           = 0.05
      ozone_data             = .false.
      co2_data               = .false.
      n2o_data               = .false.
      co_data                = .false.
      ch4_data               = .false.
      so2_data               = .false.
      addsolar               = .false.
      rayleigh_single_scatt  = .true.
      do_nlte_correction     = .false.
      solar_sea_brdf_model   = 2
      ir_sea_emis_model      = 2
      use_sfc_snow_frac      = .false.
      add_aerosl             = .false.
      aerosl_type            = 1
      add_clouds             = .true.
      ice_scheme             = 1
      use_icede              = .false.
      idg_scheme             = 2
      user_aer_opt_param     = .false.
      user_cld_opt_param     = .false.
      grid_box_avg_cloud     = .true.
      cldstr_threshold       = -1.0
      cldstr_simple          = .false.
      cldstr_low_cloud_top   = 750.0
      ir_scatt_model         = 2
      vis_scatt_model        = 1
      dom_nstreams           = 8
      dom_accuracy           = 0.0
      dom_opdep_threshold    = 0.0
      addpc                  = .false.
      npcscores              = -1
      addradrec              = .false.
      ipcreg                 = 1
      use_htfrtc             = .false.
      htfrtc_n_pc            = -1
      htfrtc_simple_cloud    = .false.
      htfrtc_overcast        = .false.
   /

| 

.. container::


   +------------------------+--------------------+----------------------------------------------------------------------+
   | Item                   | Type               | Description                                                          |
   +========================+====================+======================================================================+
   | rttov_sensor_db_file   | character(len=512) | The location of the RTTOV sensor database. The format for the        |
   |                        |                    | database is a comma-separated file. The columns of the database are  |
   |                        |                    | the DART observation-kind, the platform/satellite/sensor ID, the     |
   |                        |                    | observation type, the coefficient file, and a comma-separated list   |
   |                        |                    | of RTTOV channels to use for this observation type.                  |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | first_lvl_is_sfc       | logical            | Whether the first level of the model represents the surface (true)   |
   |                        |                    | or the top of the atmosphere (false).                                |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | mw_clear_sky_only      | logical            | If microwave calculations should be "clear-sky" only (although       |
   |                        |                    | cloud-liquid water absorption/emission is considered; see the RTTOV  |
   |                        |                    | user guide).                                                         |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | interp_mode            | integer            | The interpolation mode (see the RTTOV user guide).                   |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | do_checkinput          | logical            | Whether to check the input for reasonableness (see the RTTOV user    |
   |                        |                    | guide).                                                              |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | apply_reg_limits       | logical            | Whether to clamp the atmospheric values to the RTTOV bounds (see the |
   |                        |                    | RTTOV user guide).                                                   |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | verbose                | logical            | Whether to output lots of additional output (see the RTTOV user      |
   |                        |                    | guide).                                                              |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | fix_hgpl               | logical            | Whether the surface pressure represents the surface or the 2 meter   |
   |                        |                    | value (see the RTTOV user guide).                                    |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | do_lambertian          | logical            | Whether to include the effects of surface specularity (see the RTTOV |
   |                        |                    | user guide).                                                         |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | lambertian_fixed_angle | logical            | Whether to include a fixed angle for the lambertian effect (see the  |
   |                        |                    | RTTOV user guide).                                                   |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | rad_down_lin_tau       | logical            | Whether to use the linear-in-tau approximation (see the RTTOV user   |
   |                        |                    | guide).                                                              |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_q2m                | logical            | Whether to use 2m humidity information (see the RTTOV user guide).   |
   |                        |                    | If true, the QTY_2M_SPECIFIC_HUMIDITY will be requested from the     |
   |                        |                    | model.                                                               |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_q2m                | logical            | Whether to use 2m humidity information (see the RTTOV user guide).   |
   |                        |                    | If true, the QTY_2M_SPECIFIC_HUMIDITY will be requested from the     |
   |                        |                    | model.                                                               |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_uv10m              | logical            | Whether to use 10m wind speed information (see the RTTOV user        |
   |                        |                    | guide). If true, the QTY_10M_U_WIND_COMPONENT and                    |
   |                        |                    | QTY_10M_V_WIND_COMPONENTS will be requested from the model.          |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_wfetch             | logical            | Whether to use wind fetch information (see the RTTOV user guide). If |
   |                        |                    | true, the QTY_WIND_FETCH will be requested from the model.           |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_water_type         | logical            | Whether to use water-type information (0 = fresh, 1 = ocean; see the |
   |                        |                    | RTTOV user guide). If true, the QTY_WATER_TYPE will be requested     |
   |                        |                    | from the model.                                                      |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | addrefrac              | logical            | Whether to enable atmospheric refraction (see the RTTOV user guide). |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | plane_parallel         | logical            | Whether to treat the atmosphere as plane parallel (see the RTTOV     |
   |                        |                    | user guide).                                                         |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_salinity           | logical            | Whether to use salinity (see the RTTOV user guide). If true, the     |
   |                        |                    | QTY_SALINITY will be requested from the model.                       |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | apply_band_correction  | logical            | Whether to apply band correction from the coefficient field for      |
   |                        |                    | microwave data (see the RTTOV user guide).                           |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | cfrac_data             | logical            | Whether to use the cloud fraction from 0 to 1 (see the RTTOV user    |
   |                        |                    | guide). If true, the QTY_CLOUD_FRACTION will be requested from the   |
   |                        |                    | model.                                                               |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | clw_data               | logical            | Whether to use cloud-liquid water data (see the RTTOV user guide).   |
   |                        |                    | If true, the QTY_CLOUDWATER_MIXING_RATIO will be requested from the  |
   |                        |                    | model.                                                               |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | rain_data              | logical            | Whether to use precipitating water data (see the RTTOV user guide).  |
   |                        |                    | If true, the QTY_RAINWATER_MIXING_RATIO will be requested from the   |
   |                        |                    | model.                                                               |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | ciw_data               | logical            | Whether to use non-precipiting ice information (see the RTTOV user   |
   |                        |                    | guide). If true, the QTY_ICE_MIXING_RATIO will be requested from the |
   |                        |                    | model.                                                               |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | snow_data              | logical            | Whether to use precipitating fluffy ice (see the RTTOV user guide).  |
   |                        |                    | If true, the QTY_SNOW_MIXING_RATIO will be requested from the model. |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | graupel_data           | logical            | Whether to use precipting small, hard ice (see the RTTOV user        |
   |                        |                    | guide). If true, the QTY_GRAUPEL_MIXING_RATIO will be requested from |
   |                        |                    | the model.                                                           |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | hail_data              | logical            | Whether to use precipitating large, hard ice (see the RTTOV user     |
   |                        |                    | guide). If true, the QTY_HAIL_MIXING_RATIO will be requested from    |
   |                        |                    | the model.                                                           |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | w_data                 | logical            | Whether to use vertical velocity information. This will be used to   |
   |                        |                    | crudely classify if a cloud is cumulus or stratiform for the purpose |
   |                        |                    | of visible/infrared calculations. If true, the QTY_VERTICAL_VELOCITY |
   |                        |                    | will be requested from the model.                                    |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | clw_scheme             | integer            | The clw_scheme to use (see the RTTOV user guide).                    |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | clw_cloud_top          | real(r8)           | Lower hPa limit for clw calculations (see the RTTOV user guide).     |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | fastem_version         | integer            | Which FASTEM version to use (see the RTTOV user guide).              |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | supply_foam_fraction   | logical            | Whether to use sea-surface foam fraction (see the RTTOV user guide). |
   |                        |                    | If true, the QTY_FOAM_FRAC will be requested from the model.         |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_totalice           | logical            | Whether to use totalice instead of precip/non-precip ice for         |
   |                        |                    | microwave (see the RTTOV user guide).                                |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_zeeman             | logical            | Whether to use the Zeeman effect (see the RTTOV user guide). If      |
   |                        |                    | true, the magnetic field and cosine of bk will be used from the      |
   |                        |                    | observation metadata.                                                |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | cc_threshold           | real(r8)           | Cloud-fraction value to treat as clear-sky (see the RTTOV user       |
   |                        |                    | guide).                                                              |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | ozone_data             | logical            | Whether to use ozone (O3) profiles (see the RTTOV user guide). If    |
   |                        |                    | true, the QTY_O3 will be requested from the model.                   |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | co2_data               | logical            | Whether to use carbon dioxide (CO2) profiles (see the RTTOV user     |
   |                        |                    | guide). If true, the QTY_CO2 will be requested from the model.       |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | n2o_data               | logical            | Whether to use nitrous oxide (N2O) profiles (see the RTTOV user      |
   |                        |                    | guide). If true, the QTY_N2O will be requested from the model.       |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | co_data                | logical            | Whether to use carbon monoxide (CO) profiles (see the RTTOV user     |
   |                        |                    | guide). If true, the QTY_CO will be requested from the model.        |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | ch4_data               | logical            | Whether to use methane (CH4) profiles (see the RTTOV user guide). If |
   |                        |                    | true, the QTY_CH4 will be requested from the model.                  |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | so2_data               | logical            | Whether to use sulfur dioxide (SO2) (see the RTTOV user guide). If   |
   |                        |                    | true, the QTY_SO2 will be requested from the model.                  |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | addsolar               | logical            | Whether to use solar angles (see the RTTOV user guide). If true, the |
   |                        |                    | sun_ze and sun_az from the observation metadata will be used for     |
   |                        |                    | visible/infrared.                                                    |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | rayleigh_single_scatt  | logical            | Whether to use only single scattering for Rayleigh scattering for    |
   |                        |                    | visible calculations (see the RTTOV user guide).                     |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | do_nlte_correction     | logical            | Whether to include non-LTE bias correction for HI-RES sounder (see   |
   |                        |                    | the RTTOV user guide).                                               |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | solar_sea_brdf_model   | integer            | The solar sea BRDF model to use (see the RTTOV user guide).          |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | ir_sea_emis_model      | logical            | The infrared sea emissivity model to use (see the RTTOV user guide). |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_sfc_snow_frac      | logical            | Whether to use the surface snow fraction (see the RTTOV user guide). |
   |                        |                    | If true, the QTY_SNOWCOVER_FRAC will be requested from the model.    |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | add_aerosl             | logical            | Whether to use aerosols (see the RTTOV user guide).                  |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | aerosl_type            | integer            | Whether to use OPAC or CAMS aerosols (see the RTTOV user guide).     |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | add_clouds             | logical            | Whether to enable cloud scattering for visible/infrared (see the     |
   |                        |                    | RTTOV user guide).                                                   |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | ice_scheme             | integer            | The ice scheme to use (see the RTTOV user guide).                    |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_icede              | logical            | Whether to use the ice effective diameter for visible/infrared (see  |
   |                        |                    | the RTTOV user guide). If true, the QTY_CLOUD_ICE_DE will be         |
   |                        |                    | requested from the model.                                            |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | idg_scheme             | integer            | The ice water effective diameter scheme to use (see the RTTOV user   |
   |                        |                    | guide).                                                              |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | user_aer_opt_param     | logical            | Whether to directly specify aerosol scattering properties (see the   |
   |                        |                    | RTTOV user guide). Not yet supported.                                |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | user_cld_opt_param     | logical            | Whether to directly specify cloud scattering properties (see the     |
   |                        |                    | RTTOV user guide). Not yet supported.                                |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | grid_box_avg_cloud     | logical            | Whether to cloud concentrations are grid box averages (see the RTTOV |
   |                        |                    | user guide).                                                         |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | cldstr_threshold       | real(r8)           | Threshold for cloud stream weights for scattering (see the RTTOV     |
   |                        |                    | user guide).                                                         |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | cldstr_simple          | logical            | Whether to use one clear and one cloudy column (see the RTTOV user   |
   |                        |                    | guide).                                                              |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | cldstr_low_cloud_top   | real(r8)           | Cloud fraction maximum in layers from the top of the atmosphere down |
   |                        |                    | to the specified hPa (see the RTTOV user guide).                     |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | ir_scatt_model         | integer            | Which infrared scattering method to use (see the RTTOV user guide).  |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | vis_scatt_model        | integer            | Which visible scattering method to use (see the RTTOV user guide).   |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | dom_nstreams           | integer            | The number of streams to use with DOM (see the RTTOV user guide).    |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | dom_accuracy           | real(r8)           | The convergence criteria for DOM (see the RTTOV user guide).         |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | dom_opdep_threshold    | real(r8)           | Ignore layers below this optical depth (see the RTTOV user guide).   |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | addpc                  | logical            | Whether to do principal component calculations (see the RTTOV user   |
   |                        |                    | guide).                                                              |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | npcscores              | integer            | Number of principal components to use for addpc (see the RTTOV user  |
   |                        |                    | guide).                                                              |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | addradrec              | logical            | Reconstruct the radiances using addpc (see the RTTOV user guide).    |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | ipcreg                 | integer            | Number of predictors to use with addpc (see the RTTOV user guide).   |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | use_htfrtc             | logical            | Whether to use HTFRTC (see the RTTOV user guide).                    |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | htfrtc_n_pc            | integer            | Number of PCs to use with HTFRTC (see the RTTOV user guide).         |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | htfrtc_simple_cloud    | logical            | Whether to use simple cloud scattering with htfrtc (see the RTTOV    |
   |                        |                    | user guide).                                                         |
   +------------------------+--------------------+----------------------------------------------------------------------+
   | htfrtc_overcast        | logical            | Whether to calculate overcast radiances with HTFRTC (see the RTTOV   |
   |                        |                    | user guide).                                                         |
   +------------------------+--------------------+----------------------------------------------------------------------+

| 

Files
-----

-  A DART observation sequence file containing Radar obs.

References
----------

-  `RTTOV user guide <https://www.nwpsaf.eu/site/software/rttov/documentation/>`__

Private components
------------------

=============================== ===============================
*use obs_def_rttov_mod, only :* initialize_module
\                               initialize_rttov_sensor_runtime
\                               initialize_rttov_sensor_runtime
=============================== ===============================

| 

.. container:: routine

   *call initialize_module()*

.. container:: indent1

   Reads the namelist, allocates space for the auxiliary data associated wtih satellite observations, initializes the
   constants used in subsequent computations (possibly altered by values in the namelist), and prints out the list of
   constants and the values in use.

| 

.. container:: routine

   *call initialize_rttov_sensor_runtime(sensor,ens_size,nlevels)*
   ::

      type(rttov_sensor_type), pointer    :: sensor
      integer,                 intent(in) :: ens_size
      integer,                 intent(in) :: nlevels

.. container:: indent1

   Initialize a RTTOV sensor runtime. A rttov_sensor_type instance contains information such as options and coefficients
   that are initialized in a "lazy" fashion only when it will be used for the first time.

   ============ ===============================================
   ``sensor``   The sensor type to be initialized
   ``ens_size`` The size of the ensemble
   ``nlevels``  The number of vertical levels in the atmosphere
   ============ ===============================================

|


Error codes and conditions
--------------------------

+---------------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+
|             Routine             |                                             Message                                            |                                                              Comment                                                              |
+=================================+================================================================================================+===================================================================================================================================+
| initialize_module               | initial allocation failed for satellite observation data                                       | Need to increase MAXrttovkey                                                                                                      |
+---------------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| initialize_rttov_sensor_runtime | Module or sensor is not initialized                                                            | Both the module and the sensor must be initialized before calling this routine.                                                   |
+---------------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| get_visir_metadata              | The key exceeds the size of the metadata arrays, or the key is not a VIS/IR type               | The number of satellite observations exceeds the array size allocated in the module. Check the input and/or increase MAXrttovkey. |
+---------------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| get_mw_metadata                 | The key exceeds the size of the metadata arrays, or the key is not a MW type                   | The number of satellite observations exceeds the array size allocated in the module. Check the input and/or increase MAXrttovkey. |
+---------------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| read_rttov_metadata             | bad value for RTTOV fields                                                                     | The format of the input obs_seq file is not consistent.                                                                           |
+---------------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| get_expected_radiance           | Could not find the platform/satellite/sensor id combination in the RTTOV sensor database file. | An unknown RTTOV instrument ID was encountered. Check the database and/or the observation metadata.                               |
+---------------------------------+------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------+
 

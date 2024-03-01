.. _obs_def_rttov_mod:

MODULE ``obs_def_rttov_mod``
============================

Overview
--------

DART RTTOV observation module, including the observation operators for the two primary 
RTTOV-observation types -- visible/infrared radiances and microwave 
radiances/brightness temperatures.

DART can be built with either RTTOV v12 *or* v13. Edit :ref:`&preprocess_nml <preprocess>` to select
the appropriate obs_def:

- obs_def_rttov_mod.f90 for v12.3
- obs_def_rttov13_mod.f90 for v13 (contributed by Lukas Kugler of the University of Vienna).  

Note the namelist options for &obs_def_rttov_nml differ for v12 and v13.

- `RTTOV v12 Namelist`_ &obs_def_rttov_nml
- `RTTOV v13 Namelist`_ &obs_def_rttov_nml

For more detail on RTTOV see the `RTTOV user guide <https://www.nwpsaf.eu/site/software/rttov/documentation/>`__.

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

For RTTOV 13 DART has a ``wfetch_value`` namelist option. This allows you to set a wfetch value
to use when ``use_wfetch = .true.`` if the model you are using does not provide QTY_WIND_FETCH.

Although a model may not have the necessary inputs by itself,
the defaults in RTTOV based on climatology can be used.
The impact on the quality of the results should be investigated.

The quanities for each observation type are defined in obs_def_rttov{13}_mod.f90, like so:

.. code::

   ! HIMAWARI_8_AHI_RADIANCE,      QTY_RADIANCE

If you want to change the quantity associated with an observation, for example, if you want
to assimilate HIMAWARI_8_AHI_RADIANCE as QTY_BRIGHTNESS_TEMPERATURE, edit the QTY
in obs_def_rttov{13}_mod.f90 and rerun quickbuild.sh.


Known issues:

-  DART does not yet provide any type of bias correction
-  Cross-channel error correlations are not yet supported. A principal component approach has been discussed. For now,
   the best bet is to use a subset of channels that are nearly independent of one another.
-  Vertical localization will need to be tuned. Turning off vertical localization may work well if you have a large
   number of ensemble members. Using the maximum peak of the weighting function or the cloud-top may be appropriate.
   There are also other potential approaches being investigated.


The namelist ``&obs_def_rttov_mod_nml`` is read from file ``input.nml``. Namelists start with an ampersand '&'
and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

RTTOV v12 Namelist
------------------

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


RTTOV v13 namelist
------------------

.. code-block:: text

  &obs_def_rttov_nml
     first_lvl_is_sfc     = .true.   ! is level 1 the surface (true) or top of atmosphere (false)?
     mw_clear_sky_only    = .false.  ! only use clear-sky for MW (plus clw emission if clw_data is true) or full RTTOV-SCATT (false)?
     interp_mode          = 1        ! Interpolation mode: Rochon on OD (1), Log-linear (2), Rochon on log-linear OD (3), Rochon on WF (4), Rochon on log-linear WF (5)
     do_checkinput        = .true.   ! check if profiles are within absolute and regression limits
     apply_reg_limits     = .false.  ! clamp to min/max values
     verbose              = .true.   ! if false, only fatal errors output 
     fix_hgpl             = .true.   ! surface elevation assigned to 2m pressure (true) or surface pressure (true)
     do_lambertian        = .false.  ! treat surface as Lambertian instead of specular? (all)
     lambertian_fixed_angle = .true. ! use fixed angle for Lambertian calculations? (all, do_lambertian only)
     rad_down_lin_tau     = .true.   ! use linear-in-tau approximation? (all)
     max_zenith_angle     = 75.      ! maximum zenith angle to accept (in degrees) (all)
     use_q2m              = .false.  ! use surface humidity? (all)
     use_uv10m            = .false.  ! use u and v 10 meters? (all, used in sea surface emissivity and BRDF models)
     use_wfetch           = .false.  ! use wind fetch (length of water wind has blown over in m)  (all, used in sea surface BRDF models)
     use_water_type       = .false.  ! use water type (0 = fresh, ocean = 1) (all, used in surface BRDF atlas and models)
     addrefrac            = .true.   ! enable atmospheric refraction (all) 
     plane_parallel       = .false.  ! treat atmosphere as strictly plane-parallel? (all)
     use_salinity         = .false.  ! use ocean salinity (in practical salinity units) (MW, FASTEM 4-6 and TESSEM2)
     cfrac_data           = .false.  ! specify cloud fraction? (VIS/IR/MW)
     clw_data             = .false.  ! specify non-precip cloud liquid water? (VIS/IR/MW)
     rain_data            = .false.  ! specify precip cloud liquid water? (VIS/IR/MW)
     ciw_data             = .false.  ! specify non-precip cloud ice? (VIS/IR)
     snow_data            = .false.  ! specify precip cloud fluffy ice? (VIS/IR/MW)
     graupel_data         = .false.  ! specify precip cloud soft-hail? (VIS/IR/MW)
     hail_data            = .false.  ! specify precip cloud hard-hail? (VIS/IR/MW)
     w_data               = .false.  ! specify vertical velocity (used for classifying clouds as cumulus versus stratus)? (VIS/IR)
     clw_scheme           = 2        ! Liebe (1) or Rosenkranz (2) or TKC (3) (MW, clear-sky only)
     clw_cloud_top        = 322.0_r8   ! lower hPa limit for clw calculations; clw at lower pressures is ignored (MW, clear-sky only)
     fastem_version       = 6        ! MW sea-surface emissivity model to use (0-6). 1-6: FASTEM version 1-6, 0: TESSEM2 (MW)
     supply_foam_fraction = .false.  ! include foam fraction in skin%foam_fraction? FASTEM only. (MW)
     use_totalice         = .false.  ! Specify totalice instead of precip/non-precip ice (MW, RTTOV-SCATT only)
     use_zeeman           = .false.  ! Simulate Zeeman effect (MW)
     cc_threshold         = 0.001_r8   ! if effective cloud fraction below this value, treat simulation as clear-sky (MW, 0-1, RTTOV-SCATT only)
     ozone_data           = .false.  ! specify ozone profiles? (VIS/IR)
     co2_data             = .false.  ! specify CO2 profiles? (VIS/IR)
     n2o_data             = .false.  ! specify N2O profiles? (VIS/IR)
     co_data              = .false.  ! specify CO profiles? (VIS/IR)
     ch4_data             = .false.  ! specify CH4 profiles? (VIS/IR)
     so2_data             = .false.  ! specify SO2 profiles? (VIS/IR)
     addsolar             = .false.  ! include solar calculations (VIS/IR)
     rayleigh_single_scatt = .true.  ! if false, disable Rayleigh (VIS, addsolar only)
     do_nlte_correction   = .false.  ! if true include non-LTE bias correction for hires sounders (VIS/IR)
     solar_sea_brdf_model = 2        ! JONSWAP (1) or Elfouhaily (2) (VIS)
     ir_sea_emis_model    = 2        ! ISEM (1) or IREMIS (2) (IR)
     use_sfc_snow_frac    = .false.  ! use sfc snow cover (0-1) (IR, used in emis atlas)
     add_aerosl           = .false.  ! enable aerosol scattering (VIS/IR)
     aerosl_type          = 1        ! OPAC (1) or CAMS (2) (VIS/IR, add_aerosl only)
     add_clouds           = .true.   ! enable cloud scattering (VIS/IR)
     ice_scheme           = 1        ! SSEC (1) or Baran 2014 (2) or Baran 2018 (3) (VIS/IR, add_clouds only)
     use_icede            = .false.  ! use ice effective diameter (IR, add_clouds, ice_scheme = 1) 
     idg_scheme           = 2        ! Ou and Liou (1), Wyser (2), Boudala (3), McFarquar (2003) (VIS/IR, add_clouds only, ice_scheme = 1)
     user_aer_opt_param   = .false.  ! specify aerosol scattering properties (VIS/IR, add_clouds only)
     user_cld_opt_param   = .false.  ! specify cloud scattering properties (VIS/IR, add_clouds only)
     grid_box_avg_cloud   = .true.   ! cloud concentrations are grid box averages. False = concentrations for cloudy layer only. (VIS/IR, add_clouds and not user_cld_opt_param only)
     cldcol_threshold     = -1.0_r8    ! threshold for cloud stream weights for scattering (VIS/IR, add_clouds only)
     cloud_overlap        = 1        ! default: 1 (max/random overlap)
     cc_low_cloud_top     = 750.0_r8   ! cloud fraction maximum in layers from ToA down to specified hPa (VIS/IR, cloud_overlap only)
     ir_scatt_model       = 2        ! DOM (1) or Chou-scaling (2) (IR, add_clouds or add_aerosl only)
     vis_scatt_model      = 1        ! DOM (1), single scat (2), or MFASIS (3) (VIS, addsolar and add_clouds or add_aerosl only)
     dom_nstreams         = 8        ! number of streams to use with DOM (VIS/IR, add_clouds or add_aerosl and DOM model only, must be >= 2 and even)
     dom_accuracy         = 0.0_r8     ! convergence criteria for DOM (VIS/IR, add_clouds or addaerosol and DOM model only)
     dom_opdep_threshold  = 0.0_r8     ! DOM ignores layers below this optical depth (VIS/IR, add_clouds or addaerosol and DOM model only)
     addpc                = .false.  ! do principal component calculations? (VIS/IR)
     npcscores            = -1       ! number of PC scores to use (VIS/IR, addpc only)
     addradrec            = .false.  ! reconstruct the radiances (VIS/IR, addpc only)
     ipcreg               = 1        ! number of predictors, see Table 29 of user guide (VIS/IR, addpc only)
     use_htfrtc           = .false.  ! use HTFRTC of Havemann 2018  
     htfrtc_n_pc          = -1       ! number of PCs to use (HTFRTC only, max 300)
     htfrtc_simple_cloud  = .false.  ! use simple-cloud scattering (HTFRTC only)
     htfrtc_overcast      = .false.  ! calculate overcast radiances (HTFRTC only)
     wfetc_value          = 100000.0_r8 ! Real wfetc Wind fetch (m) (length of water over which the wind has blown, typical
                                                              ! value 100000m for open ocean). Used if wfetc not provided by model.
  /

References
----------

-  `RTTOV user guide <https://www.nwpsaf.eu/site/software/rttov/documentation/>`__



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
 

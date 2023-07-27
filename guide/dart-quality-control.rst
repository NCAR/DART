DART quality control field
==========================

DART has a quality control (QC) field in the *obs_seq.final* file
to report on the status of the assimilation of the variable. The most
common reason for exploring the DART QC value is to help determine
if the observation was assimilated (or evaluated) - or if the observation
was rejected or ...

To learn more about how to intepret the QC field as well as
other values in an observation sequence file,
see :doc:`detailed-structure-obs-seq`.
The ‘DART QC’ field is usually the second of the 2 “quality control” copies.

A list of all the DART QC values can be found in the QC table in
:doc:`../assimilation_code/modules/assimilation/quality_control_mod`.


-  If the DART QC values are 4, the forward operators have failed. Look at the
   *model_interpolate()* routine in your model_mod.f90 file, or the forward
   operator code in *observations/forward_operators/obs_def_xxx_mod.f90* for
   your observation type. A successful forward operator must return a valid
   obs_val and an *istatus = 0*. If the forward operator code returns different
   istatus values for different error types, you can set
   *&filter_nml::output_forward_op_errors = .true.* and rerun *filter* to see
   exactly what error istatus codes are being set. See
   :doc:`../assimilation_code/modules/assimilation/filter_mod` for more
   information on how to use the ‘output_forward_op_errors’ option. Negative
   istatus values are reserved for the system, *istatus = 0* is success, and any
   positive value indicates a failed forward operator. The code is free to use
   different positive values to signal different types of errors.

-  If the DART QC values are 5, those observation types were intentionally
   ignored because they were not listed in the &obs_kind_nml namelist, in the
   ‘assimilate_these_obs_types’ stringlist.

-  If the DART QC values are 6, the data quality control that came with the
   original observation data indicates this is a bad quality observation and it
   was skipped for this reason.

-  If the DART QC values are 7, the observation value is too far away from the
   ensemble mean. Set *&quality_control_nml::outlier_threshold = -1* to ignore this for
   now and rerun. In general, this is not the optimal strategy as the number of
   observations inconsistent with the ensemble is a very powerful indicator of
   the success or failure of the assimilation.

-  If the DART QC values are 8, it was not possible to convert the observation
   to the required vertical coordinate system.

If the prior and posterior values in the ``obs_seq.final`` are not -888888.0 but
are identical, your obs are being assimilated but are having no impact.

The most common reasons assimilated obs have no impact on the model state
include:

-  **Zero spread in ensemble members**
   Your initial ensemble members must have different values for each state item.
   If all members have identical values, the observations cannot make a change.
   To diagnose this condition, look at the prior ensemble spread. This is either
   in ``preassim.nc`` or ``preassim_sd.nc``, depending on your model. If all the
   values are 0, this is your problem. One way to generate an ensemble with some
   spread is to set *&filter_nml::perturb_from_single_instance = .false.,*
   (which will still require a single filter initial condition file) but then
   the *filter* code will add random gaussian perturbations to each state vector
   item to generate an initial ensemble with spread. The magnitude of the
   gaussian noise added is controlled by the
   *&filter_nml::perturbation_amplitude*. It is also possible to write your own
   perturbation routine in your ``model_mod.f90`` code.
-  **Cutoff value too small**
   If the localization radius is too small, the observation may not be ‘close
   enough’ to the model grid to be able to impact the model. Check the
   localization radius (*&assim_tools_nml::cutoff*). Set it to a very large
   number (e.g. 100000) and rerun. If there is now an impact, the cutoff was
   restricting the items in the state vector so your obs had no impact before.
   Cutoff values are dependent on the location type being used. It is specified
   in radians for the threed_sphere locations module (what most large models
   use), or in simple distance (along a unit circle) if using a low order model
   (lorenz, ikeda, etc).
-  **Obs error values too large (less likely)**
   If the observation error is very large, it will have no impact on the model
   state. This is less likely a cause than other possibilities.
-  **No correlation (unlikely)**
   If there is no correlation between the distribution of the forward
   observation values and the state vector values, the increments will be very
   tiny. However there are generally still tiny increments applied, so this is
   also a low likelyhood case.
-  **Errors in forward operator location computations, or get_close_obs()**
   If there is an error in the ``model_mod.f90`` code in either
   *get_state_meta_data()*, *model_interpolate()*, or the vertical conversion
   code in *get_close_obs()*, it is possible for the forward operators to appear
   to be working correctly, but the distances computed for the separation
   between the obs and the state vector values can be incorrect. The most
   frequent problem is that the wrong locations are being passed back from
   *get_state_meta_data()*. This can result in the increments being applied in
   the wrong locations or not at all. This is usually one of the things to test
   carefully when developing a new model interface, and usually why we recommend
   starting with a single observation at a known location.
-  **Incorrect vertical conversion**
   If the model is using 3d coordinates and needs the capability to convert
   between pressure, height, and/or model level, the conversion may be
   incorrect. The state vector locations can appear to be too high or too low to
   be impacted by an observation. Some models have a height limit built into
   their model_mod code to avoid trying to assimilate observations at the model
   top. The observations cannot make meaningful changes to the model state there
   and trying to assimilate them can lead to problems with the inflation. If the
   code in the model_mod is excluding observations incorrectly, or you are
   testing with observations at the model top, this can result in no impact on
   the model state.

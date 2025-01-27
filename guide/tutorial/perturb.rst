Quick Note on Using the Perturbed Observation EnKF
==================================================

DART supports sorting of observation increments before computing state increments.

This option should not be used with any filter except the perturbed observation EnKF.
It should always be used with the EnKF to improve filter results.

The perturbed observation ensemble Kalman filter can be selected by using the enkf_qceff_table.csv.

Use the enkf table in the namelist and set sort_obs_inc = .true. in the assim_tools_nml section.
Run filter and looks at the results. 

Then change back to the bnrhf_qceff_table.csv and set sort_obs_inc = .false.
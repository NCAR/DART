Runtime control via input.nml in the work directory
===================================================

The file input.nml in the work directory allows the value of many DART parameters to be set at runtime. 

The input.nml contains a section for each DART module and can include controls for multiple main programs 
like perfect_model_obs and filter. 

For example, input.nml contains the namelist variables for the DART module assim_tools_nml:

.. code-block:: text

	&assim_tools_nml
	   cutoff                          = 0.2,
	   sort_obs_inc                    = .false.,
	   spread_restoration              = .false.,
	   sampling_error_correction       = .false.,
	   adaptive_localization_threshold = -1,
	   output_localization_diagnostics = .false.,
	   localization_diagnostics_file   = 'localization_diagnostics',
	   print_every_nth_obs             = 0,
	   rectangular_quadrature          = .true.,
	   gaussian_likelihood_tails       = .false.,
	/

The following input.nml variables allow exploration of DART filter features that were discussed in the 
earlier DART_LAB tutorial sections.

1. In the assim_tools_nml, cutoff: Sets the halfwidth of a Gaspari-Cohn localization for observation impact.

2. In the filter_nml, ens_size: Sets the ensemble size for the assimilation.

See how changing the ensemble size and cutoff impact the assimilation results.

- Change localization by setting cutoff to 0.4. Run filter and use matlab tools to see impact.
- Change ens_size to 80 (keep 0.4 cutoff). Check the impacts.


Inflation is controlled by a block of entries in the filter_nml section of input.nml:

.. code-block:: text

	inf_flavor                  = 0,         0,
	inf_initial_from_restart    = .false.,   .false.,
	inf_sd_initial_from_restart = .false.,   .false.,
	inf_deterministic           = .true.,    .true.,
	inf_initial                 = 1.0,       1.0,
	inf_lower_bound             = 0.0,       1.0,
	inf_upper_bound             = 1000000.0, 1000000.0,
	inf_damping                 = 0.9,       1.0, 
	inf_sd_initial              = 0.6,       0.0,
	inf_sd_lower_bound          = 0.6,       0.0,
	inf_sd_max_change           = 1.05,      1.05,

The entries in the first column of numbers control prior inflation, that is applied after the model 
advance and before the assimilation. This is what was available in DART_LAB.

DART also supports posterior inflation, controlled by the second column of numbers, which applies 
inflation after the assimilation but before the next model advance. Prior and posterior inflation 
is also supported.

The first row in the inflation namelist controls is inf_flavor. This value controls the algorithmic 
variant of inflation applied. The following values are currently supported:

.. list-table:: Inflation Flavor
   :header-rows: 1

   * - Flavor
     - Description
   * - 0
     - No inflation
   * - 1
     - Time varying adaptive inflation with a single value for all variables at a given time
   * - 2
     - Spatially- and temporally-varying using a Gaussian prior for the inflation value
   * - 5
     - Spatially- and temporally-varying with inverse gamma prior for the inflation value


Try changing the prior inf_flavor (first column) from 0 to 5 (keep all other namelist settings).

Then try changing the ens_size back to 20. 

The combination of inflation, localization, and ensemble size controls the assimilation quality.

Other inflation controls that were discussed in DART_LAB include:

.. list-table:: Inflation Controls
   :header-rows: 1

   * - Control
     - Description
   * - inf_lower_bound
     - Inflation is not allowed to be smaller than this value.
   * - inf_upper_bound
     - Inflation is not allowed to exceed this value.
   * - inf_damping
     - The inflation is damped towards 1 by this factor at each assimilation time.
   * - inf_sd_initial
     - The inflation standard deviation initial value.
   * - inf_sd_lower_bound
     - Inflation lower bound cannot be smaller than this.
   * - inf_sd_max_change
     - Fractional change in inflation standard deviation cannot exceed this at a given assimilation time.


The values for the prior inflation (column 1) set in the default input.nml in the lorenz_96 work directory 
are a good choice for many applications.  
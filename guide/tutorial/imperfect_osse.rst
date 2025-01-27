Imperfect Model OSSEs
=====================

Parameters of the forecast model are set in the model_nml

.. code-block:: fortran

	&model_nml
		model_size = 40,
		forcing = 8.00,
		delta_t  = 0.05,
		time_step_days = 0,
		time_step_seconds = 3600,
	/

The forcing, F, of the Lorenz_96 model is set to 8 by default.
Changing this value when running perfect_model_obs or filter leads to different forcing.
Running perfect_model_obs with one value (like 8) and the filter with another explores 
assimilation with an imperfect model as in the matlab gui.

Set the model forcing to 10

.. code-block:: fortran

	&model_nml
		model_size        = 40,
		forcing           = 10.00,
		delta_t           = 0.05,
		time_step_days    = 0,
		time_step_seconds = 3600,
	/

Run filter and see how the results have been degraded.

Change the forcing back to 8 and run filter again before moving ahead.
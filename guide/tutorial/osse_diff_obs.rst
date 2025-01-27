OSSEs with Different Observations
=================================

The observation sequence used to generate observations for an OSSE is set by:

.. code-block:: text

	&perfect_model_obs_nml
	 	...
		obs_seq_in_file_name	= “obs_seq.in”

To do an OSSE with a different observing network:

- Change this file name to one of the others on the previous slide,
- Run perfect_model_obs to generate synthetic obs, then run filter.
- Try this with obs_seq_identity_1_20.in.
- Use plot_ens_time_series to look at results at different points.

Change back to the default observation sequence file, obs_seq.in, before moving ahead:

.. code-block:: text

	&perfect_model_obs_nml
	    ...
		obs_seq_in_file_name	= “obs_seq.in”

Run perfect_model_obs to recreate the original set of synthetic observations.

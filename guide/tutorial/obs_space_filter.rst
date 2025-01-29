Controlling the observation space filter algorithm
==================================================

Some DART namelist entries are the names of files that contain more detailed run-time control information. 

The most widely used example is found in the algorithm_info module namelist:

.. code-block:: text

	&algorithm_info_nml
	   qceff_table_filename = 'eakf_qceff_table.csv'
	/

The default namelist entry in Lorenz_96/work is the eakf_qceff_table.csv

In this case, DART uses an Ensemble Adjustment Kalman Filter in observation space and does no 
transforms before regressing increments. This is equivalent to using a normal distribution.
.. _quantile tracer:

QCEFF: Examples with the Lorenz 96 Tracer Model
===============================================


The Quantile-Conserving Ensemble Filter Framework (QCEFF) tools are available in DART
as of version v11.
The DART development team (dart@ucar.edu) would be happy to hear about your experiences and is
anxious to build scientific collaborations using these new capabilities.

To get started, make sure that you are on the quantile_methods branch of DART: 

.. code-block:: text

   git checkout quantile_methods

Build the DART executables for the Lorenz 96 tracer advection model:

.. code-block:: text

    cd DART/models/lorenz_96_tracer_advection/work
    ./quickbuild.sh nompi


The new quantile options are set using a :ref:`qceff table <QCEFF>` given as a namelist
option ``qceff_table_filename`` to &algorithm_info_nml. The examples below show how to change the quantile options
using various QCEFF tables. You can find the .csv files for these four examples in the directory
``DART/models/lorenz_96_tracer_advection/work``


.. list-table::
   :header-rows: 1 
   :widths: 15 60 25

   * - example
     - description
     - .csv filename 
   * - Example A 
     - bounded normal rank histogram
     - all_bnrhf_qceff_table.csv
   * - Example B
     - Ensemble Adjustment Kalman filters
     - all_eakf_qceff_table.csv 
   * - Example C
     - EAKF for state and bounded normal rank histogram filter and priors for tracer concentration and source
     - state_eakf_tracer_bnrhf_qceff_table.csv
   * - Example D
     - Negative tracers bounded above
     - neg_qceff_table.csv


You can view .csv files with a text editor, or spreadsheet tool such as Google Sheets,
or Microsoft Excel.

Example A
----------

Assimilating observations of state (wind) and tracer concentration using
a rank histogram observation space filter and rank histogram probit transforms for
state variable updates. This example includes adaptive inflation.

The default model configuration has a single tracer source at gridpoint 1 along with
small uniform tracer sinks that lead to areas where the true tracer concentration is
usually 0. This is a particularly tough test for ensemble methods.

#. Edit input.nml to set the qceff_table_filename to 'all_bnrhf_qceff_table.csv' 

   .. code-block:: text

      &algorithm_info_nml
         qceff_table_filename = 'all_bnrhf_qceff_table.csv'
       

#. Create a set_def.out file using create_obs_sequence,

   ``./create_obs_sequence < create_obs_sequence_input``

#. Create an obs_sequence.in file using create_fixed_network_seq

      ``./create_fixed_network_seq``

   .. code:: text

      Select the default input filename <return>,
      Create a regularly repeating sequence by entering "1",
      Enter "1000" for the number of observation times,
      Enter "0 0" for the initial time,
      Enter "0 10800" for the period,
      Select the default output filename, <return>

#. Spin-up a model initial condition by running perfect_model_obs

      ``./perfect_model_obs``

#. Generate a spun-up true time series,

      ``cp perfect_output.nc perfect_input.nc``


   Edit input.nml to set read_input_state_from_file to .true.

   .. code:: text
     
      &perfect_model_obs_nml
        read_input_state_from_file = .true.,


   Run ``./perfect_model_obs`` again.

#. Run a filter assimilation,

      ``./filter``

#. Examine the output with your favorite tool(s) (e.g. plot_ens_time_series.m). Looking at the analysis ensemble 
   for the tracer_concentration variables with indices near the source (location 1)
   and far downstream from the source (location 35) is interesting.
   Near the source, the true concentration and the ensemble estimates are all non-zero while far from the source
   there are times when the true concentration and many ensemble members are zero. For further detail
   see Anderson et al. (2023). [1]_
   Note that the source estimation capabilities of the model and filters are not being tested here.


Example B 
---------

Using Ensemble Adjustment Kalman filters.


#. Edit input.nml to set the qceff_table_filename to 'all_eakf_qceff_table.csv'

   .. code-block:: text

      &algorithm_info_nml
         qceff_table_filename = 'all_eakf_qceff_table.csv'
       

#. Run the filter 

      ``./filter``

Example C 
---------

Using Ensemble Adjustment Kalman filter for state, but bounded normal rank histogram filter and priors for tracer concentration and source.


#. Edit input.nml to set the qceff_table_filename to state_eakf_tracer_bnrhf_qceff_table.csv

   .. code-block:: text

      &algorithm_info_nml
         qceff_table_filename = 'state_eakf_tracer_bnrhf_qceff_table.csv'
       

#. Run the filter 

     ``./filter``

Example D 
----------

Testing the bounded above option. Normally tracers are bounded below, but there are other quantities that may be bounded
above. There are distinct numerical challenges in implementing the quantile algorithms
for quantities that are bounded above, so flipping the sign of the tracers is a good
test. 

#. Edit input.nml to set the qceff_table_filename to neg_qceff_table.csv

   .. code-block:: text

      &algorithm_info_nml
         qceff_table_filename = 'neg_qceff_table.csv'
      

#. Edit input.nml, to change the entry positive_tracer to .false. and read_input_state_from_file back to .false. 

   
   .. code-block:: text

      &model_nml
          positive_tracer          = .false.,

      &perfect_model_obs_nml
          read_input_state_from_file = .false.,


#. Repeat steps 3-6 from Test A.

References
----------

.. [1] Anderson, J. L., Riedel, C., Wieringa, M., Ishraque, F., Smith, M., Kershaw, H.
       2023: A Quantile-Conserving
       Ensemble Filter Framework. Part III: Data Assimilation for Mixed Distributions
       with Application to a Low-Order Tracer Advection Model. *Monthly Weather Review*
       `[Manuscript submitted for publication] <../_static/papers/QCEFF_3_submitted.pdf>`_

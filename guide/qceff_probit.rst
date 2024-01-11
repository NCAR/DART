.. _QCEFF:

Quantile-Conserving Ensemble Filter Framework
==============================================

The Quantile-Conserving Ensemble Filter Framework (QCEFF) tools are available in DART
as of version v11. 
The DART development team (dart@ucar.edu) would be happy to hear about your experiences 
and is anxious to build scientific collaborations using these new capabilities.

The QCEFF options are set using a :ref:`qceff table <qceff table>` given as a namelist option to &algorithm_info_nml.

   .. code-block:: text

      &algorithm_info_nml
         qceff_table_filename = 'qceff_table.csv'


.. _QCEFF options:

QCEFF options
--------------

QCEFF options are per quantity. For a given quantity, you specify the following 
options as columns of the qceff_table:

* Observation error information

   Provides information about boundedness constraints that control the likelihood
   distribution associated with an observed variable when using perfect_model_obs
   to generate noisy observations.
   
     * bounded_below (default .false.) 
     * bounded_above (default .false.)
     * lower_bound   
     * upper_bound


* Probit distribution information 

   Used in the computation of the probit transform.
   The values needed are the bounds and the distribution type.
   These can be different for all three cases (inflation, state, and extended_state priors)
   
     * distribution (one of :ref:`Distributions`)
     * bounded_below (default .false.)
     * bounded_above (default .false.)
     * lower_bound    (default -888888)
     * upper_bound    (default -888888)


* Observation increment information

     * filter_kind (one of :ref:`Filter kinds`)
     * bounded_below (default .false.)
     * bounded_above (default .false.)
     * lower_bound    (default -888888)
     * upper_bound    (default -888888)



.. _qceff table:

Creating a qceff table
-----------------------

The table has two headers, row 1 and 2.
The first row is the version number.  The second row describes the :ref:`QCEFF options` corresponding to each column of the table. 
These two headers must be present in your qceff table. 
The :ref:`qcf trunc table` below shows the table headers, 
and an example quantity QTY_TRACER_CONCENTRATION for the first 5 columns of the table. 
There is a complete table with all 25 columns in `Google Sheets <https://docs.google.com/spreadsheets/d/1CRGHWc7boQt81pw2pDxEFY6WPyQeCh64OwPyoVMqijE/edit?usp=sharing>`_. You can copy and edit this table as needed.

To add a quantity, add a row to the table.
For any quantity not listed in the table, the :ref:`Default values` values will be used for all 25 options. 
You only have to add rows for quantities that use non-default values for any of the input options.
Ensure that there are no empty rows in between the quantities listed in the spreadsheet.
Save your spreadsheet as a .csv file. 

To run filter or perfect_model_obs, put the .csv file in the directory where you are running.
Edit input.nml to set the algorithm_info_nml option qceff_table_filename, for example:


   .. code-block:: text

      &algorithm_info_nml
         qceff_table_filename = 'qceff_table.csv'


.. _qcf trunc table:

.. list-table:: truncated table 
   :header-rows: 2

   * - QCF table version: 1
     - 
     -  
     -  
     - 
   * - QTY
     - bounded_below
     - bounded_above
     - lower_bound
     - upper_bound
   * - QTY_TRACER_CONCENTRATION
     - .true.
     - .false.
     - 0
     - -888888


.. _Filter kinds:

Available filter kinds
-----------------------

   * EAKF (default)
   * ENKF
   * UNBOUNDED_RHF
   * GAMMA_FILTER
   * BOUNDED_NORMAL_RHF

.. _Distributions:

Available distributions
------------------------

  * NORMAL_DISTRIBUTION (default)
  * BOUNDED_NORMAL_RH_DISTRIBUTION
  * GAMMA_DISTRIBUTION 
  * BETA_DISTRIBUTION
  * LOG_NORMAL_DISTRIBUTION
  * UNIFORM_DISTRIBUTION



.. _Default values:

Default values
---------------

If a quantity is not in the qceff table, the following default values
are used:

  * filter_kind (default EAKF)
  * dist_type (default NORMAL_DISTRIBUTION)
  * bounded_below  (default .false.)
  * bounded_above   (default .false.)
  * lower_bound    (default -888888)
  * upper_bound    (default -888888)

.. note::

   -888888 is a missing value in DART.


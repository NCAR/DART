MARBL_column
============

**MARBL** stands for the Marine Biogeochemistry Library; it's a modular biogeochemical modeling suite for next-generatioon models. 
It simulates marine ecosystem dynamics and the coupled cycles of carbon, nitrogen, phosphorus, iron, silicon, and oxygen. 
It is a component of the Community Earth System Model (`CESM <https://www.cesm.ucar.edu/>`_) and is often coupled to 
other physical ocean models such as the Modular Ocean Model (`MOM6 <https://mom6.readthedocs.io/en/main/>`_). 
                         
For a detailed description of the model, the reader is refered to Long et al., 2021 [1]_.
                         
The MARBL source code can be downloaded from https://github.com/marbl-ecosys/MARBL. You may also want to read 
MARBL's documentation `here <https://marbl-ecosys.github.io>`_.

This MARBL-DART interface was developed by `Robin Armstrong <https://github.com/robin-armstrong>`_. Thanks Robin! 

Overview 
--------
**MARBL_column** provides an interface between DART and a 1D (column) configuration of MARBL. 
It's designed to be used at locations where in-situ data is abundant such as 
`Bermuda Atlantic Time-series Study (BATS) station <https://bats.bios.asu.edu/>`_, 
`Weather Station Mike <https://projects.met.no/iaoos/en/en-testitem1/work-packages/wp3-process-experiments/task-3.2-towards-a-modern-weather-station-mike/indexccb4.html>`_, 
`Hawaii Ocean Time-series (HOT) <https://hahana.soest.hawaii.edu/hot/>`_ ... to name a few. 
 
The code is designed to perform 3 kinds of data assimilation (DA) experiments: 
                         
#. **State Estimation:** where the prognostic state variables of MARBL such as nitrate concerntration are updated.
   To achive this, you'll need to set  
   
   ``estimate_params = .false.`` within ``&model_nml`` in the namelist file ``input.nml`` 
                         
#. **State and Parameters Estimation:** where both the state and a set of model parameters are updated. 
   MARBL has a long list of uncertain model parameters that can be constrained alongside the state. 
   This usually improves the prediction skill of the model and alleviates some of its biases. 
   To achieve this DA exercise, you'll need to set
   
   ``estimate_params = .true.`` 
   
   The combined DART state will be of the form :math:`Z_k = \left[ \mathbf{x}_k, \boldsymbol{\theta} \right]^T`
   where :math:`Z_k` is the joint state, :math:`\mathbf{x}` and parameters, :math:`\boldsymbol{\theta}` 
   vector at time :math:`t_k`      
                    
#. **Parameters Estimation only:** where only the parameters are constrained using the data. DART
   will still need to read in the state to construct ensemble covariances and compute innovations. 
   The only difference is that the ensemble increments are only regressed onto the unknown parameters.
   To achive this goal, you'll need to set ``estimate_params = .true.`` and turn the update status for 
   all state variable to ``NO_COPY_BACK``
   
   This ensures that the updated state will not be written back to the restart file. The ``NO_COPY_BACK`` 
   option is added as the 5th entry in the state table (after the variable name, its associated quantity 
   and its physical bounds) within ``&model_nml``. 

Namelist
--------
The ``&model_nml`` variables and their default values are listed here:

.. code-block:: fortran 

  &model_nml
     state_template_file   = 'MOM.res.nc', 
     param_template_file   = 'marbl_params.nc',
     ocean_geometry        = 'ocean_geometry.nc',
     station_location      = -64, 31
     time_step_days        = 1,
     time_step_seconds     = 0,
     model_state_variables = 'NO3      ', 'QTY_NITRATE_CONCENTRATION     ', '0.0', 'NA', 'UPDATE      ',
                             'SiO3     ', 'QTY_DISSOLVED_INORGANIC_SIO3  ', '0.0', 'NA', 'UPDATE      ',
                             'PO4      ', 'QTY_PHOSPHATE_CONCENTRATION   ', '0.0', 'NA', 'UPDATE      ',
                             'Fe       ', 'QTY_DISSOLVED_INORGANIC_IRON  ', '0.0', 'NA', 'UPDATE      ',
                             'DIC      ', 'QTY_DISSOLVED_INORGANIC_CARBON', '0.0', 'NA', 'UPDATE      ',
                             'O2       ', 'QTY_DISSOLVED_OXYGEN          ', '0.0', 'NA', 'UPDATE      ',
                             'DOC      ', 'QTY_DISSOLVED_ORGANIC_CARBON  ', '0.0', 'NA', 'UPDATE      ',
                             'DON      ', 'QTY_DISSOLVED_ORGANIC_NITROGEN', '0.0', 'NA', 'UPDATE      ',
                             'DOP      ', 'QTY_DISSOLVED_ORGANIC_P       ', '0.0', 'NA', 'UPDATE      ',
                             'ALK      ', 'QTY_ALKALINITY                ', '0.0', 'NA', 'UPDATE      ',
                             'microzooC', 'QTY_MICROZOOPLANKTON_CARBON   ', '0.0', 'NA', 'UPDATE      ',
                             'mesozooC ', 'QTY_MESOZOOPLANKTON_CARBON    ', '0.0', 'NA', 'UPDATE      ',
                             'h        ', 'QTY_LAYER_THICKNESS           ', '0.0', 'NA', 'NO_COPY_BACK'
     estimate_params       = .true.
     model_parameters      = 'autotroph_settings(1)%kDOP ', 'QTY_BGC_PARAM', '0.0', 'NA', 'UPDATE'
     /

This namelist provides control over the kind of DA experiment as described abvove. 

+-------------------------------------+--------------------+------------------------------------------------------------+
| Item                                | Type               | Description                                                |
+=====================================+====================+============================================================+
| ``state_template_file``             | character(len=256) | MARBL restart file including MARBL's prognostic variables  |
|                                     |                    | and other grid information such as the ocean layers.       |
|                                     |                    | The state variables read from this file are listed in      |
|                                     |                    | in the ``model_state_variables``                           |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``param_template_file``             | character(len=256) | Template file corresponding to the BGC parameters that we  |
|                                     |                    | intend to estimate. The parameters read from this file are |
|                                     |                    | listed in the ``model_parameters``.                        |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``ocean_geometry``                  | character(len=256) | The ocean geometry file is used to read ``basin_depth``    |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``station_location``                | real(2)            | Longitude and latitude of the ocean column location.       |
+-------------------------------------+--------------------+------------------------------------------------------------+                               
| ``time_step_days``                  | integer            | The number of days to advance the model for each           | 
|                                     |                    | assimilation.                                              |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``time_step_seconds``               | integer            | In addition to ``time_step_days``, the number              |
|                                     |                    | of seconds to advance the model for each assimilation.     |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``model_state_variables``           | character(:,5)     | Strings that associate MARBL variables with a DART         |
|                                     |                    | quantity. They also describe their physical bounds and     |
|                                     |                    | whether or not to write the updated values to the restart  |
|                                     |                    | files. These variables will be read from the MARBL restart |
|                                     |                    | file and modified by the assimilation. Some (perhaps all)  |
|                                     |                    | will be used by the forward observation operators. If the  |
|                                     |                    | 5th column is ``UPDATE``, the output files will have the   |
|                                     |                    | modified (assimilated,posterior) values. If the 5th        |
|                                     |                    | column is ``NO_COPY_BACK``, that variable will not be      |
|                                     |                    | written to the restart files. **The DART diagnostic files  |
|                                     |                    | will always have the (modified) posterior values.**        |
|                                     |                    | Diagnostic variables that are useful for the calculation   |
|                                     |                    | of the forward observation operator but have no impact on  |
|                                     |                    | the forecast trajectory of the model could have a value of |
|                                     |                    | ``NO_COPY_BACK``. The 3rd and 4th column list the minimum  |
|                                     |                    | and maximum allowed values for the updated variables.      |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``estimate_params``                 | logical            | A switch to turn on/off parameter estimation.              |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``model_parameters``                | character(:,5)     | Similar to ``model_state_variables``, this is a list of    |
|                                     |                    | parameters that will take part in the DART state and       |
|                                     |                    | would possibly get updated.                                |
+-------------------------------------+--------------------+------------------------------------------------------------+


References
----------
.. [1] Long, Matthew C., J. Keith Moore, Keith Lindsay, Michael Levy, Scott C. Doney, 
       Jessica Y. Luo, Kristen M. Krumhardt, Robert T. Letscher, Maxwell Grover, and Zephyr T. Sylvester. 
       "Simulations with the marine biogeochemistry library (MARBL)." 
       Journal of Advances in Modeling Earth Systems 13, no. 12 (2021): e2021MS002647.

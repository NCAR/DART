null_model
==========

Overview
--------

DART interface module for the 'null_model'. This model provides very simple
models for evaluating filtering algorithms. It can provide simple linear growth
around a fixed point, a random draw from a Gaussian, or combinations of the two.
Namelist controls can set the width of the Gaussian and change both the model
advance method and the expected observation interpolation method.

The 18 public interfaces are standardized for all DART compliant models. These
interfaces allow DART to advance the model, get the model state and metadata
describing this state, find state variables that are close to a given location,
and do spatial interpolation for model state variables.

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: fortran

   &model_nml
      model_size           = 2,
      delta_t              = 0.05,
      time_step_days       = 0,
      time_step_seconds    = 3600  
      noise_amplitude      = 0.0_r8
      advance_method       = 'simple'
      interpolation_method = 'standard'
   /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------------+---------------+--------------------------------------+
| Item                 | Type          | Description                          |
+======================+===============+======================================+
| model_size           | integer       | Model size.                          |
+----------------------+---------------+--------------------------------------+
| delta_t              | real(r8)      | Internal model timestep parameter.   |
+----------------------+---------------+--------------------------------------+
| time_step_days       | integer       | Minimum model advance time in days.  |
+----------------------+---------------+--------------------------------------+
| time_step_seconds    | integer       | Minimum model advance time in        |
|                      |               | seconds.                             |
+----------------------+---------------+--------------------------------------+
| noise_amplitude      | real(r8)      | If greater than 0.0 sets the         |
|                      |               | standard deviation of the added      |
|                      |               | Gaussian noise during the model      |
|                      |               | advance.                             |
+----------------------+---------------+--------------------------------------+
| advance_method       | character(64) | Controls the model advance method.   |
|                      |               | The default is 'simple'              |
|                      |               | timestepping. A 4-step Runga Kutta   |
|                      |               | method can be selected with the      |
|                      |               | string 'rk'.                         |
+----------------------+---------------+--------------------------------------+
| interpolation_method | character(64) | Controls how the expected value of   |
|                      |               | an observation is computed. The      |
|                      |               | default is 'standard' which uses a   |
|                      |               | linear interpolation between the two |
|                      |               | surrounding model points. Other      |
|                      |               | options include 'square' which       |
|                      |               | returns the square of the computed   |
|                      |               | value, 'opposite_side' which adds on |
|                      |               | a value from the opposite side of    |
|                      |               | the cyclical domain, and 'average'   |
|                      |               | which averages 15 points to get the  |
|                      |               | expected value. Model size should be |
|                      |               | > 15 to use the last option.         |
+----------------------+---------------+--------------------------------------+

Files
-----

+-----------------------------+-----------------------------------------------+
| filename                    | purpose                                       |
+=============================+===============================================+
| input.nml                   | to read the model_mod namelist                |
+-----------------------------+-----------------------------------------------+
| preassim.nc                 | the time-history of the model state before    |
|                             | assimilation                                  |
+-----------------------------+-----------------------------------------------+
| analysis.nc                 | the time-history of the model state after     |
|                             | assimilation                                  |
+-----------------------------+-----------------------------------------------+
| dart_log.out [default name] | the run-time diagnostic output                |
+-----------------------------+-----------------------------------------------+
| dart_log.nml [default name] | the record of all the namelists actually USED |
|                             | - contains the default values                 |
+-----------------------------+-----------------------------------------------+

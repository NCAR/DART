Controlling which files are output by filter
============================================

DART provides you with fine-grained control over how and when files are output.
You can instruct DART whether or not to output files after each stage in an
assimilation cycle. Since most experiments are run for more than one
assimilation cycle, you can also instruct DART to aggregate all of the output
for a specific stage into a single file.

These options are controlled by three settings in the ``filter_nml`` namelist
in ``input.nml``:

- ``stages_to_write`` specifies the stages during an assimilation cycle during
  which state files may be output. The possible stages are 
  ``'input'``, ``'forecast'``, ``'preassim'``,
  ``'postassim'``, ``'analysis'`` and ``'output'``. The input strings are
  case-insensitive, but the corresponding output files are always lowercase.
- ``single_file_in`` specifies how input state files are structured. If 
  ``.true.`` the state of all ensemble members is expected to be read from
  single file. If  ``.false.`` the state of each ensemble member expected to 
  be read from its own file.
- ``single_file_out`` specifies how output state files are structured. If 
  ``.true.`` the state of all ensemble members is output to a single file. If 
  ``.false.`` the state of each ensemble members is output to its own file.

.. caution::

  ``single_file_out`` only refers to the output **for a particular stage**.
  So even if you set ``single_file_out = .true.``, you can get *several* output
  files - one per stage. If you set ``single_file_out = .false.`` filter will
  output a deluge of files. Be careful about what stages you choose to write.

Two common assimilation workflows
---------------------------------

There are many ways to configure your data assimilation workflows. However, the
following two workflows are sensible for small models and large models,
respectively.

Small models
~~~~~~~~~~~~

For models that read and write small state files and complete their numerical
integrations relatively quickly, it makes sense to configure filter to:

1. complete multiple assimilation cycles
2. read from and write to a single output file for all ensemble members

This workflow requires setting ``single_file_in = .true.`` and
``single_file_out = .true.``.

When ``filter`` is used for a long assimilation experiment, setting
``single_file_out = .true.`` will consolidate all the information for a
particular stage into a single file that contains all the ensemble members,
the mean, spread, inflation, etc.

This results in far fewer files, and each file may contain multiple timesteps
to encompass the entirety of the experiment. Take note: since a single task
must write each file, this setting engenders some computational overhead.

Large models
~~~~~~~~~~~~

For models that read and write large state files and complete their numerical
integrations relatively slowly, it make sense to configure filter to:

1. complete a single assimilation cycle at a time
2. read from and write to a seperate output file for each ensemble member

This workflow requires setting ``single_file_out = .false.`` and makes sense 
for large models or in cases where it is beneficial to run different number of
MPI tasks for the model advances and the assimilation. In this case, there can
be a substantial computational efficiency to have each ensemble member write
its information to a separate file, and each file can be written simultaneously
by different tasks. The tradeoff (at the moment) is that each of the files can
only have a single timestep in them. Consequently, some files are redundant and
should not be output.

Output and diagnostic files produced by filter
----------------------------------------------

In the case when ``single_file_out = .false.``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

======================== ========== =================================================================
from *perfect_model_obs*            
======================== ========== =================================================================
``obs_seq.out``                     the synthetic observations at some predefined times and locations
``perfect_output.nc``    1 timestep  a netCDF file containing the model trajectory - the true state
======================== ========== =================================================================

There are some namelist settings that control what files are output. Depending on the settings for
*input.nml&filter_nml:stages_to_write* and others …

+--------------------------------------+---------------------+---------------------+
| from *filter*                        |                     |                     |
+======================================+=====================+=====================+
| ``forecast_member_####.nc``          | 1 timestep          | the ensemble        |
|                                      |                     | forecast, each      |
|                                      |                     | ensemble member is  |
|                                      |                     | a separate file     |
+--------------------------------------+---------------------+---------------------+
| ``forecast_[mean,sd].nc``            | 1 timestep          | the mean and        |
|                                      |                     | standard deviation  |
|                                      |                     | (spread) of the     |
|                                      |                     | ensemble forecast   |
+--------------------------------------+---------------------+---------------------+
| ``forecast_priorinf_[mean,sd].nc``   | 1 timestep          | the prior inflation |
|                                      |                     | information before  |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``forecast_postinf_[mean,sd].nc``    | 1 timestep          | the posterior       |
|                                      |                     | inflation           |
|                                      |                     | information before  |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``preassim_member_####.nc``          | 1 timestep          | the model states    |
|                                      |                     | after any prior     |
|                                      |                     | inflation but       |
|                                      |                     | before assimilation |
+--------------------------------------+---------------------+---------------------+
| ``preassim_[mean,sd].nc``            | 1 timestep          | the mean and        |
|                                      |                     | standard deviation  |
|                                      |                     | (spread) of the     |
|                                      |                     | ensemble after any  |
|                                      |                     | prior inflation but |
|                                      |                     | before assimilation |
+--------------------------------------+---------------------+---------------------+
| ``preassim_priorinf_[mean,sd].nc``   | 1 timestep          | the prior inflation |
|                                      |                     | information before  |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``preassim_postinf_[mean,sd].nc``    | 1 timestep          | the posterior       |
|                                      |                     | inflation           |
|                                      |                     | information before  |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``postassim_member_####.nc``         | 1 timestep          | the model states    |
|                                      |                     | after assimilation  |
|                                      |                     | but before          |
|                                      |                     | posterior inflation |
+--------------------------------------+---------------------+---------------------+
| ``postassim_[mean,sd].nc``           | 1 timestep          | the mean and        |
|                                      |                     | standard deviation  |
|                                      |                     | (spread) of the     |
|                                      |                     | ensemble after      |
|                                      |                     | assimilation but    |
|                                      |                     | before posterior    |
|                                      |                     | inflation           |
+--------------------------------------+---------------------+---------------------+
| ``postassim_priorinf_[mean,sd].nc``  | 1 timestep          | the (new) prior     |
|                                      |                     | inflation           |
|                                      |                     | information after   |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``postassim_postinf_[mean,sd].nc``   | 1 timestep          | the (new) posterior |
|                                      |                     | inflation           |
|                                      |                     | information after   |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``analysis_member_####.nc``          | 1 timestep          | the model states    |
|                                      |                     | after assimilation  |
|                                      |                     | and after any       |
|                                      |                     | posterior inflation |
+--------------------------------------+---------------------+---------------------+
| ``analysis_[mean,sd].nc``            | 1 timestep          | the mean and        |
|                                      |                     | standard deviation  |
|                                      |                     | (spread) of the     |
|                                      |                     | ensemble after      |
|                                      |                     | assimilation and    |
|                                      |                     | after posterior     |
|                                      |                     | inflation           |
+--------------------------------------+---------------------+---------------------+
| ``analysis_priorinf_[mean,sd].nc``   | 1 timestep          | the (new) prior     |
|                                      |                     | inflation           |
|                                      |                     | information after   |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``analysis_postinf_[mean,sd].nc``    | 1 timestep          | the (new) posterior |
|                                      |                     | inflation           |
|                                      |                     | information after   |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``output_[mean,sd].nc``              | 1 timestep          | the mean and spread |
|                                      |                     | of the posterior    |
|                                      |                     | ensemble            |
+--------------------------------------+---------------------+---------------------+
| ``output_priorinf_[mean,sd].nc``     | 1 timestep          | the (new) prior     |
|                                      |                     | inflation           |
|                                      |                     | information after   |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``output_priorinf_[mean,sd].nc``     | 1 timestep          | the (new) posterior |
|                                      |                     | inflation           |
|                                      |                     | information after   |
|                                      |                     | assimilation        |
+--------------------------------------+---------------------+---------------------+
| ``obs_seq.final``                    |                     | the model estimates |
|                                      |                     | of the observations |
|                                      |                     | (an integral part   |
|                                      |                     | of the data         |
|                                      |                     | assimilation        |
|                                      |                     | process)            |
+--------------------------------------+---------------------+---------------------+

+---------------------+-----------------------------------+
| from both           |                                   |
+=====================+===================================+
| ``dart_log.out``    | the ‘important’ run-time output   |
|                     | (each run of *filter* appends to  |
|                     | this file; remove it or start at  |
|                     | the bottom to see the latest      |
|                     | values)                           |
+---------------------+-----------------------------------+
| ``dart_log.nml``    | the input parameters used for an  |
|                     | experiment                        |
+---------------------+-----------------------------------+

In the case when ``single_file_out = .true.``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All the information for each stage is contained in a single file that *may* have multiple timesteps.

======================== =========== =================================================================
from *perfect_model_obs*             
======================== =========== =================================================================
``obs_seq.out``                      the synthetic observations at some predefined times and locations
``perfect_output.nc``    N timesteps a netCDF file containing the model trajectory - the true state
======================== =========== =================================================================

There are some namelist settings that control what files are output. Depending on the settings for ``input.nml
&filter_nml:stages_to_write`` and others.

+------------------------+--------------------+--------------------+
| from *filter*          |                    |                    |
+========================+====================+====================+
| ``filter_input.nc``    | 1 timestep         | The starting       |
|                        |                    | condition of the   |
|                        |                    | experiment. All    |
|                        |                    | ensemble members,  |
|                        |                    | [optionally] the   |
|                        |                    | input mean and     |
|                        |                    | standard deviation |
|                        |                    | (spread),          |
|                        |                    | [optionally] the   |
|                        |                    | prior inflation    |
|                        |                    | values,            |
|                        |                    | [optionally] the   |
|                        |                    | posterior          |
|                        |                    | inflation values   |
+------------------------+--------------------+--------------------+
| ``forecast.nc``        | N timesteps        | The ensemble       |
|                        |                    | forecast. All      |
|                        |                    | ensemble members,  |
|                        |                    | the mean and       |
|                        |                    | standard deviation |
|                        |                    | (spread), the      |
|                        |                    | prior inflation    |
|                        |                    | values, the        |
|                        |                    | posterior          |
|                        |                    | inflation values   |
+------------------------+--------------------+--------------------+
| ``preassim.nc``        | N timesteps        | After any prior    |
|                        |                    | inflation but      |
|                        |                    | before             |
|                        |                    | assimilation. All  |
|                        |                    | ensemble members,  |
|                        |                    | the mean and       |
|                        |                    | standard deviation |
|                        |                    | (spread) of the    |
|                        |                    | ensemble, the      |
|                        |                    | prior inflation    |
|                        |                    | values, the        |
|                        |                    | posterior          |
|                        |                    | inflation values   |
+------------------------+--------------------+--------------------+
| ``postassim.nc``       | N timesteps        | After assimilation |
|                        |                    | but before         |
|                        |                    | posterior          |
|                        |                    | inflation. All     |
|                        |                    | ensemble members,  |
|                        |                    | the mean and       |
|                        |                    | standard deviation |
|                        |                    | (spread) of the    |
|                        |                    | ensemble, the      |
|                        |                    | (new) prior        |
|                        |                    | inflation values,  |
|                        |                    | the (new)          |
|                        |                    | posterior          |
|                        |                    | inflation values   |
+------------------------+--------------------+--------------------+
| ``analysis.nc``        | N timesteps        | After assimilation |
|                        |                    | and after any      |
|                        |                    | posterior          |
|                        |                    | inflation. All     |
|                        |                    | ensemble members,  |
|                        |                    | the mean and       |
|                        |                    | standard deviation |
|                        |                    | (spread) of the    |
|                        |                    | ensemble, the      |
|                        |                    | (new) prior        |
|                        |                    | inflation values,  |
|                        |                    | the (new)          |
|                        |                    | posterior          |
|                        |                    | inflation values   |
+------------------------+--------------------+--------------------+
| ``filter_output.nc``   | 1 timestep         | After assimilation |
|                        |                    | and after any      |
|                        |                    | posterior          |
|                        |                    | inflation. All     |
|                        |                    | ensemble members,  |
|                        |                    | the mean and       |
|                        |                    | standard deviation |
|                        |                    | (spread) of the    |
|                        |                    | ensemble, the      |
|                        |                    | (new) prior        |
|                        |                    | inflation values,  |
|                        |                    | the (new)          |
|                        |                    | posterior          |
|                        |                    | inflation values   |
+------------------------+--------------------+--------------------+
| ``obs_seq.final``      |                    | the model          |
|                        |                    | estimates of the   |
|                        |                    | observations (an   |
|                        |                    | integral part of   |
|                        |                    | the data           |
|                        |                    | assimilation       |
|                        |                    | process)           |
+------------------------+--------------------+--------------------+

+-----------------------+---------------------+
| from both             |                     |
+=======================+=====================+
| ``dart_log.out``      | the ‘important’     |
|                       | run-time output     |
|                       | (each run of        |
|                       | *filter* appends to |
|                       | this file; remove   |
|                       | it or start at the  |
|                       | bottom to see the   |
|                       | latest values)      |
+-----------------------+---------------------+
| ``dart_log.nml``      | the input           |
|                       | parameters used for |
|                       | an experiment       |
+-----------------------+---------------------+

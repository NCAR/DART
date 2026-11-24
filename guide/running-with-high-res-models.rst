.. _high-res-tips:

Tips for Running with High Resolution Models
============================================

One of the more challenging aspects of an ensemble Data Assimilation (DA) system 
is the need to manage large amounts of memory to store ensembles of the model state.

Most modern large-scale models run in parallel on multi-processor computer systems, 
distributing data across multiple memory nodes to support finer spatial grids, 
smaller time steps, and longer simulation periods. Common strategies include using 
shared memory within individual nodes and the Message Passing Interface (MPI) 
on distributed-memory systems.

Ensemble DA exacerbates this memory problem by requiring multiple copies,
often 20-100x, of the model data to do the assimilation.

DART uses the Message Passing Interface (MPI) to distribute ensembles of 
state and extended state data across processors and nodes in a distributed-memory system. 
For models with relatively small memory demands, DART can be compiled and run as a serial program. 
However, when compiled with MPI, it can scale to tens of thousands of processors 
and handle memory requirements ranging from gigabytes to terabytes. This allows DART to run with 
extremely high resolution models and large ensembles. 

There is an inherent trade-off between memory usage and internode communication.
DART provides multiple runtime strategies to balance this trade-off: 
users can choose to reduce memory usage per MPI task at the cost of increased communication
overhead, or to use more memory per node in order to minimize data transfer between MPI tasks.

It may be necessary to adjust these runtime settings when working with high resolution or large models, 
as the memory requirements may exceed the memory available on a single node. This document will detail 
the different namelist options that control the distribution of data across MPI tasks, as well as 
provide recommendations for requesting resources on multi-node systems to enable working with larger 
model states.

For quick reference, here is a collection of the recommended namelist settings for running with 
high resolution models, which are discussed in more detail below. Ensure that the following 
namelist settings are set appropriately and fit your system's resources and assimilation needs.

.. code-block:: text

    &ensemble_manager_nml
       layout                      = 2
       tasks_per_node              = 128 ! Note this should match your mpirun settings; 128 is the recommended 
                                         ! value when using all available MPI tasks per node on Derecho
    /
                                         
    &filter_nml
       distributed_state = .true.,
       stages_to_write          = 'output',
       single_file_out          = .false.,
       output_members           = .true.,
       output_mean              = .true.,
       output_sd                = .true.,
       num_output_state_members = 0,         ! or a small subset of the ensemble
       write_all_stages_at_end  = .false.,
       compute_posterior        = .false.,   ! switch to .true. if you want posterior inflation or need O-A diagnostics
       num_output_obs_members   = 0          ! or a small subset of the ensemble
    /

    &obs_sequence_nml
       write_binary_obs_sequence = .false.   ! switch to .true. if obs_seq.final is too large to write in ASCII
    /

    &assim_tools_nml
       distribute_mean = .false.
    /

    &state_vector_io_nml
       buffer_state_io          = .true.
    /


Distributed State
~~~~~~~~~~~~~~~~~

A state vector is read into memory in a single MPI task, and then distributed across all the
MPI tasks. By default, DART runs in distributed mode, ``distributed_state = .true.``,
where each MPI task has all the ensemble members for subset of the state vector. The entire 
filter run keeps this data distribution and the data is only collected back into individual 
ensembles at output time. For the phases in DART that require access to the whole state vector, 
DART takes care of retrieving any state values needed from other MPI tasks. 

Recommended setting for ``distributed_state``:

.. code-block:: text

    &filter_nml
       distributed_state = .true.,
    /

.. note::

   For debugging purposes, you may want to run in "transpose mode". If all data values for a single 
   ensemble member are needed, data transposes are used to collect the state vector data onto the first N 
   tasks, where N is the ensemble size.

   This is done by setting the namelist option ``distributed_state = .false.`` in the 
   filter_nml namelist. In general, this is not recommended as it requires the full state vector 
   to fit into a single task’s memory, has less vectorization, and does not scale beyond N tasks. 

To minimize per node memory and IO contention, we recommend changing some parameters in the ensemble_manager_nml
namelist. 

The first is the ``layout`` parameter, which controls the mapping between logical MPI tasks 
and physical MPI tasks. We recommend setting this to a value of 2, which tells DART to distribute the 
logical MPI tasks across the physical MPI tasks in a round-robin fashion. This means that each node gets 
a task in turn before assigning a second task to each node, until all tasks are assigned. This helps to 
balance memory usage across nodes and can improve performance.

If the ``layout`` parameter is not set in the ensemble_manager_nml namelist, DART will use the default 
value of 1, which is a one-to-one mapping between logical and physical MPI tasks (PE = MPI task).

Once you have set the ``layout`` parameter to 2, you must also set the ``tasks_per_node`` parameter so that 
the number of logical MPI tasks per node matches the number of physical MPI tasks per node you are using 
on your machine.

If ``tasks_per_node`` is not set in the ensemble_manager_nml namelist, DART will use the default 
value of 1 MPI task per node. The job will still run if this number does not match the real “tasks_per_node” at the hardware level, 
but it may run out of memory if the mismatch causes multiple high-memory-use tasks to be run on the 
same node.

Recommended setting for a round-robin layout in ensemble_manager_nml:

.. code-block:: text

    &ensemble_manager_nml
       layout                      = 2
       tasks_per_node              = 128 ! Note this should match your mpirun settings; 128 is the recommended 
                                         ! value when using all available MPI tasks per node on Derecho


IO - Reading and Writing of the Model State and Observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first N logical MPI tasks read and write the N state vectors from disk, where N is the ensemble size. 
For very large state vectors, you may not be able to read the entire state vector into memory on a single 
node. In this case you can set ``buffer_state_io = .true.`` in the state_vector_io_nml namelist.
This will read/write the state vector in chunks, buffering the data in memory before distributing it across 
the MPI tasks.

    .. code-block:: text

        &state_vector_io_nml
           buffer_state_io          = .true.,
        /

.. warning::

   Buffering state IO has more MPI communication overhead than the default setting of ``buffer_state_io = .false.``,
   so it is only recommended if you have a very large state vector that does not fit into memory on a single node.


DART gives you fine-grained control over which output files are written and when. You can choose which stages of
the assimilation cycle writes state files, how many ensemble members go into each file, and how much
information is written to the observation-space diagnostic file, ``obs_seq.final``.

For small models, writing everything costs little, so these settings aren't critical. For large models,
writing unnecessary output can take a significant share of the runtime, use more memory, and fill up disk
space. It helps to know which files the assimilation actually needs and which are only for diagnostics:

- **Needed:** the input (prior) ensemble, which is the set of model restart files that filter reads, and the
  updated (posterior) ensemble, which filter writes back to the model restart files in the ``'output'`` stage.
  The model needs these to continue to the next cycle.
- **Diagnostic only:** the ``'forecast'``, ``'preassim'``, ``'postassim'`` and ``'analysis'`` stages, the
  posterior forward operators, and the ensemble members written to ``obs_seq.final``. These do not change
  the assimilation result. However, they provide extra information and data to analyze and can be useful 
  for checking and tuning an experiment. For example, comparing ``forecast`` and ``preassim`` members shows 
  the effect of prior inflation. 

For large models, it is recommended that you write only what you need for routine cycling, and turn on the 
diagnostic output only when you need it, for example for a short test period.

These options are set in the ``&filter_nml`` namelist in ``input.nml`` unless noted otherwise.

``stages_to_write`` (default: ``'output'``)
   The stages of the assimilation cycle at which state files are written. Valid values are ``'input'``,
   ``'forecast'``, ``'preassim'``, ``'postassim'``, ``'analysis'``, ``'output'`` and ``'null'``. The values
   are case-insensitive, but the output file names are always lowercase.

   - ``'input'``: the mean and spread of the input ensemble (the members themselves are already the input files).
   - ``'forecast'``: the ensemble before prior inflation is applied.
   - ``'preassim'``: after prior inflation, before assimilation.
   - ``'postassim'``: after assimilation, before posterior inflation.
   - ``'analysis'``: after assimilation and posterior inflation.
   - ``'output'``: the updated ensemble written to the restart files for the next model advance. All
     ensemble members are always written in this stage.

   For large models, ``stages_to_write = 'output'`` is usually all that is needed. Each additional stage
   writes more files, which adds I/O time and disk use.

``output_mean`` (default: ``.true.``) and ``output_sd`` (default: ``.true.``)
   Whether the ensemble mean and ensemble standard deviation (spread) are written in
   each enabled stage. The mean and spread are each the size of a single ensemble member, so they are typcailly 
   cheap to write. Unless you are explicitly running into memory constraints, it is recommended to 
   keep these namelist settings to ``.true.``, as they are useful for diagnostics and tuning.

``output_members`` (default: ``.true.``)
   Whether the ensemble members are written in each enabled stage.

``num_output_state_members`` (default: ``0``)
   The number of ensemble members written in the ``'forecast'``, ``'preassim'``, ``'postassim'`` and
   ``'analysis'`` stages. ``output_members`` must be ``.true.``. With the default of 0, these stages write only
   the mean, spread, and inflation values. Writing members in these stages multiplies the output size, so for
   large models keep this at 0 or a small number. This setting does not affect the ``'output'`` stage.

``single_file_out`` (default: ``.false.``)
   How output state files are structured. If ``.true.``, all ensemble members for a stage are written to a single
   file. If this is set to ``.false.``, each member is written to its own file, with the first N tasks each 
   concurrently write the output file for a single ensemble member. This is the recommended setting for large models.

``write_all_stages_at_end`` (default: ``.false.``)
   If ``.false.``, state data for each stage is written as soon as it is computed. If ``.true.``, filter keeps a
   copy of every requested stage (members, mean, spread, and inflation) in memory and writes all of them in
   parallel at the end of the run. This can lower the cost of jobs where I/O time dominates, but it increases
   memory use for every stage written. It cannot be used with ``single_file_out = .true.``. For memory-limited
   runs, keep the default ``.false.``.

``compute_posterior`` (default: ``.true.``)
   Whether forward operators are computed on the posterior ensemble. If computed, these values are written to
   ``obs_seq.final``. They do not change the assimilation result. Setting this to ``.false.`` skips the 
   posterior forward operator calculation, saving time and memory, and halves the number of ensemble copies 
   written to ``obs_seq.final``. It is recommended to set this to ``.false.`` for large models, unless you explicitly 
   need the posterior forward operators for specific diagnostics or posterior inflation.

``num_output_obs_members`` (default: ``0``)
   The number of ensemble members whose forward operator values are written to ``obs_seq.final``, in addition to
   the ensemble mean and spread. Each member adds one value per observation for the prior, and another for the
   posterior if ``compute_posterior = .true.``. With many observations, writing the whole ensemble can make
   ``obs_seq.final`` extremely large. MPI task 0 also holds all of these values in memory and writes the file by
   itself, so it is recommended to write a subset of the ensemble, or use 0 if you only need the mean and spread.

``write_binary_obs_sequence`` in ``&obs_sequence_nml`` (default: ``.false.``)
   If ``.true.``, observation sequence files, including ``obs_seq.final``, are written in binary instead of
   ASCII. Binary files are smaller and faster to read and write, which helps when an ASCII ``obs_seq.final``
   becomes too large. They are not human-readable. If you read them on a machine with a different byte order,
   set ``read_binary_file_format`` in ``&obs_sequence_nml``. Recommended to keep this setting ``.false.`` 
   unless you are running into memory issues with ASCII files when writing the final observation sequence.

To summarize, these are the recommended settings for large, memory-limited models:

.. code-block:: text

    &filter_nml
       stages_to_write          = 'output',
       single_file_out          = .false.,
       output_members           = .true.,
       output_mean              = .true.,
       output_sd                = .true.,
       num_output_state_members = 0,         ! or a small subset of the ensemble
       write_all_stages_at_end  = .false.,
       compute_posterior        = .false.,   ! switch to .true. if you want posterior inflation or need O-A diagnostics
       num_output_obs_members   = 0          ! or a small subset of the ensemble
    /

    &obs_sequence_nml
       write_binary_obs_sequence = .false.    ! switch to .true. if obs_seq.final is too large in ASCII
    /

.. note::

   Turn on the diagnostic stages, ensemble members, and posterior forward operators when you need them to
   evaluate an experiment, for example for a short period while tuning inflation and localization.


Conversion of Observation or State Location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The vertical location of an observation or state element may depend on the model state. For example, 
an observation might be reported in pressure coordinates, requiring conversion to, for example, height,
based on the model state.
For localization calculations, all MPI tasks use the ensemble mean state for any vertical coordinate 
transformations that depend on the model state.

If your state is small or your system has ample memory, it is advantageous to store a full copy of the 
mean on each MPI task for use in vertical calculations. However, if the state is very large and you're 
approaching (or exceeding) the available memory per node, you can reduce memory usage by distributing the 
mean across MPI tasks.

The default for ``distribute_mean`` is ``.false.``, and it is recommended to keep it as such for most cases. 
However, if you are running into memory constraints (the memory of an entire state vector can not fit into 
a single task’s memory) or you do not need to perform vertical coordinate conversions, set this to ``.true.``.

Note that distributing the mean can significantly hurt your performance by increasing MPI communication overhead,
so only use it when necessary.

.. code-block:: text

    &assim_tools_nml
       distribute_mean = .true.,
    /


Advice for Running on Multi-node Machines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have tried all of the above recommendations and are still running into memory issues, you can use the technique of
"under-subscribing" compute nodes to reduce the number of MPI tasks that are assigned per node. This will increase the amount 
of memory that is available for each MPI task.

For example, on NSF NCAR's supercomputer Derecho, each CPU node has 256 GB of memory. If you are running a job with 
128 MPI tasks per node (utilizing the maximum CPUs for computation), each MPI task will have 2 GB of memory available. 
If you were to reduce the number of MPI tasks per node to 64, it will result in each MPI task having 4 GB of memory available.

.. note:: 

   When altering the number of MPI tasks per node, you must also remember to change the ``tasks_per_node`` 
   parameter in the ensemble_manager_nml namelist.

   .. code-block:: text

      &ensemble_manager_nml
         tasks_per_node              = 64 ! Note this should match your mpirun settings


The file system you are running on can make a significant difference in the performance of DART, especially 
when dealing with high resolution models or large state vectors. 

For example, on NSF NCAR's supercomputer Derecho, the centralized file service uses high-performance 
Spectrum Scale and Lustre shared file system technology. This allows users to run jobs on either the 
Lustre or Spectrum Scale file systems. Running on the Lustre file system (scratch) greatly reduces 
runtime due to Lustre systems having tools designed to maximize the usability and performance of both 
read/write operations and metadata operations.

DART does have strong scaling, so utlizing more nodes can greatly improve the performance of DART 
by reducing each processor's workload. However, it is important to be aware of the limits to strong 
scaling, or strong scaling breakdown. 

There is a limit to how much adding more nodes will improve performance. As you add more processors 
to a fixed problem, you eventaully hit a point of diminishing returns where communication overwhelms 
computation and your runtime will start to increase. It is important to be aware of this limit and to 
avoid using more nodes than necessary, as it can lead to diminished performance, wasted resources, 
and increased energy costs.
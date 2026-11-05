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
and handle memory requirements ranging from gigabytes to petabytes. This allows DART to run with 
extremely high resolution models and large ensembles. 

There is an inherent trade-off between memory usage and internode communication.
DART provides multiple runtime strategies to balance this trade-off: 
users can choose to reduce memory usage per MPI task at the cost of increased communication
overhead, or to use more memory per node in order to minimize data transfer between MPI tasks.

It may be necessary to adjust these runtime settings when working with high resolution or large models, 
as the memory requirements may exceed the memory available on a single node. This document will detail 
the different namelist options that control the distribution of data across MPI tasks, as well as 
provide recommendations for requesting resources on multi-node systems to enable working with larger 
model states. Ensure that the following namelist settings are set appropriately.


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
   ensemble are needed, data transposes are used to collect the state vector data onto the first N 
   tasks, where N is the ensemble size.

   This is done by setting the namelist option ``distributed_state = .false.`` in the 
   filter_nml namelist. In general, this is not recommended as it requires the full state vector 
   to fit into a single task’s memory, has less vectorization, and does not scale beyond N tasks. 

To minimize per node memory and IO contention, we recommend you assign the logical MPI tasks 
to physical MPI tasks using a round-robin layout in the ensemble_manager_nml namelist.

You must also set the ``tasks_per_node`` parameter so that the number of logical MPI 
tasks per node matches the number of physical MPI tasks per node you are using on your machine.

Recommended setting for a round-robin layout in ensemble_manager_nml:

.. code-block:: text

    &ensemble_manager_nml
       layout                      = 2
       tasks_per_node              = 128 ! Note this should match your mpirun settings


IO - reading and writing of the model state
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Conversion of observation or state location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The vertical location of an observation or state element may depend on the model state. For example, 
an observation might be reported in pressure coordinates, requiring conversion to, for example, height,
based on the model state.
For localization calculations, all MPI tasks use the ensemble mean state for any vertical coordinate 
transformations that depend on the model state.

If your state is small or your system has ample memory, it can be advantageous to store a full copy of the 
mean on each MPI task for use in vertical calculations. However, if the state is very large and you're 
approaching (or exceeding) the available memory per node, you can reduce memory usage by distributing the 
mean across MPI tasks.

The default for ``distribute_mean`` is ``.false.``  but for models with large state vectors or models that 
do not perform vertical coordinate conversions, it may be beneficial to set this to ``.true.``.

.. code-block:: text

    &assim_tools_nml
       distribute_mean = .true.,
    /

.. warning::

   Unless you are explicitly running into memory constraints, it is not recommended to change this namelist setting
   to ``distribute_mean = .true.``, as it can significantly hurt your performance by increasing MPI communication overhead.


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

For a specific assimilation (i.e. fixed problem size), there is a limit to how much adding more nodes 
will improve performance. As you add more processors to a fixed problem, you eventaully hit a point 
of diminishing returns where communication overwhelms computation and your runtime will start to 
increase. It is important to be aware of this limit and to avoid using more nodes than necessary, 
as it can lead to diminished performance, wasted resources, and increased energy costs.
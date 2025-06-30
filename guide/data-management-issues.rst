.. _data-distribution:

Data management in DART	
=======================

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
and handle memory requirements ranging from gigabytes to petabytes.

There are three phases in DART that require access to the whole state vector:

1. IO - reading and writing of the model state.
2. Forward Operator (FO) computation. 
3. Conversion of an observation or state location from one vertical coordinate system to another.

There is an inherent trade-off between memory usage and internode communication.
DART provides multiple runtime strategies to balance this trade-off: 
users can choose to reduce memory usage per MPI task at the cost of increased communication
overhead, or to use more memory per node in order to minimize data transfer between MPI tasks.
The whole state is always logically visible to all MPI tasks, but physically it may be distributed 
across multiple tasks.


IO - reading and writing of the model state
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first N logical MPI tasks read and write the N state vectors from disk, where N is the ensemble size.
A state vector is read into memory in a single MPI task, and then distributed across all the
MPI tasks.  To minimize per node memory and IO contention, we recommend you assign the logical MPI tasks 
to physical MPI tasks using a round-robin layout in the ensemble_manager_nml namelist.

Recommended setting for a round-robin layout in ensemble_manager_nml:

.. code-block:: text

    &ensemble_manager_nml
       layout                      = 2
       tasks_per_node              = 128 ! Note this should match your mpirun settings
    /

For very large state vectors, you may not be able to read the entire state vector into memory on a single 
node. In this case you can set ``buffer_state_io = .true.`` in the state_vector_io_nml namelist.
This will read the state vector in chunks, buffering the data in memory before distributing it across the MPI tasks.

.. Note:: 

    Note buffering state IO has more communication overhead than the default setting of ``buffer_state_io = .false.``,
    so is only recommended if you have a very large state vector that does not fit into memory on a single node.


    .. code-block:: text

        &state_vector_io_nml
           buffer_state_io          = .true.,
        /


Forward Operator (FO) computation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The majority of calculation in an assimilation is done across the ensemble, for example,
means, variances, increments, inflation, etc. and so the ideal data layout for assimilation is to have 
all ensemble members for a given element of the state vector on the same MPI task. 
However, :ref:`forward operator <FO>` calculations may need data from any part of the state
vector, and so the ideal data layout for the forward operator is the entire state vector on 
the same MPI task.

By default, DART runs in distributed mode, ``distributed_state = .true.`` (ideal for assimilation)
where each MPI task has all the ensemble members for subset of the state vector.
The forward operator calculation is vectorized across the whole ensemble, and DART takes 
care of retrieving any state values needed for the forward operator calculation from other MPI
tasks. 

You may want to run in "transpose mode" where the first N MPI tasks perform the forward operator
calculation, where N is the ensemble size. This is done by setting the namelist option
``distributed_state = .false.`` in the filter_nml namelist. In general, this is not recommended
as it requires more memory per task, has less vectorization, and does not scale beyond N tasks. 
However, it may be useful
for debugging.

Recommended setting for ``distributed_state``:

.. code-block:: text

    &filter_nml
       distributed_state = .true.,
    /


Conversion of observation or state location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The vertical location of an observation or state element may depend on the model state. For example, 
an observation might be reported in pressure coordinates, requiring conversion to, for example, height,
based on the model state.
For localization calculations, all MPI tasks use the ensemble mean state for any vertical coordinate 
transformations that depend on the model state.

If your state is small or your system has ample memory, it can be advantageous to store a full copy of the 
mean on each MPI task for use in vertical calculations. However, if the state is very large and you're 
approaching (or exceeding) the available memory per node, it may be more efficient to distribute the mean 
across MPI tasks to reduce memory usage

The default for ``distribute_mean`` is ``.false.``  but for models with large state vectors or models that 
do not perform vertical coordinate conversions, it may be beneficial to set this to ``.true.``.

.. code-block:: text

    &assim_tools_nml
       distribute_mean = .false.,
    /





Data management in DART	
=======================

One of the more challenging aspects of an ensemble Data Assimilation (DA) system 
is the need to manage large amounts of memory to store ensembles of the model state.

Most contemporary large models run in parallel on multi-processor computer systems 
and distribute the
data across multiple memory nodes to support finer grids, smaller timesteps and
longer modeling time periods.  Common computer science strategies include using
shared memory on individual nodes and using the Message Passing Interface (MPI) 
libraries on distributed memory nodes.

Ensemble DA exacerbates this memory problem by requiring multiple copies,
often 20-100x, of the model data to do the assimilation.

DART uses the MPI libraries to distribute ensembles of model state data across
distributed memory nodes.  For models with small amounts of data the code can be 
compiled and run as a serial program but when compiled with MPI 
it can scale up to 10,000s of nodes using Giga to Petabytes of memory.

Memory usage and internode communication time are mutually incompatible 
items to minimize.  DART has different strategies that can be selected
at runtime to use less memory per node at the cost of more time spent in communication
of data between nodes, or use more memory per node and minimize communication time.

The following descriptions detail the different phases of the main assimilation
program in DART, called ``filter``, and what options exist for memory layout
and management.


Ensembles of data
-----------------

State data
~~~~~~~~~~

* N ensemble members times X items in the state vector, always resident.
* 6 additional copies of X items for inflation, ensemble mean & sd, etc.

Observations
~~~~~~~~~~~~

Allocated and deallocated if looping over multiple assimilation windows
within a single run of filter.

* Only observations within the current assimilation window, O
* O observations times N ensemble members for the Forward Operator (FO) results
* O observations times N ensemble members for the QC results

Delayed writing option
~~~~~~~~~~~~~~~~~~~~~~

If selected in the namelist, up to P phases (input, forecast, preassim, postassim,
analysis, output) of the state data are stored in memory and written out at the end 
of filter.


Filter run phases
-----------------

FO computation, prior and posterior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run-time options include allocating spaces for two layouts and transposing
between them, or running distributed in 'all copies' mode.


Assimilation
~~~~~~~~~~~~

Distributed FO and QC observation ensembles

Runtime option to either replicate the model state ensemble mean on each MPI task
or run with that ensemble fully distributed.


Ensemble memory usage and layout
--------------------------------

Transposable
~~~~~~~~~~~~

Data is distributed over T MPI tasks but during the program execution
the data is communicated between tasks to alternate between two different
data layouts.

Allocations are needed for two different 2D arrays: 

* N ensemble members times (X items/T tasks)
* X items times (N ensemble members/T tasks)

Distributed
~~~~~~~~~~~

Data is distributed over T MPI tasks but only a single data array
is used:

* N ensemble members times (X items/T tasks)

Replicated
~~~~~~~~~~

The same data array is replicated on each MPI task:

* X items per task

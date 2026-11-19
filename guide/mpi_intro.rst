Message Passing Interface
=========================

Introduction
------------

DART programs can be compiled using the *Message Passing Interface (MPI)*.
MPI is both a library and run-time system that enables multiple copies of a
single program to run in parallel, exchange data, and combine to solve a
problem more quickly.

DART does **NOT** require MPI to run. However, for larger models with large state vectors and
large numbers of observations, the data assimilation step will run much faster
in parallel, which requires MPI to be installed and used. However, if multiple
ensembles of your model fit comfortably (in time and memory space) on a single
processor, you need read no further about MPI.

MPI is a `standard <https://www.mpi-forum.org/>`_; there are many implementations of MPI. If you
have a large single-vendor system it probably comes with an MPI library by
default. For a Linux cluster there are generally more variations in what might
be installed; examples include `OpenMPI <https://www.open-mpi.org/>`_ and 
`MVAPICH2 <http://mvapich.cse.ohio-state.edu/>`_.


An "MPI program" makes calls to an MPI library, and needs to be compiled with MPI include files and libraries.
Generally the MPI installation includes a shell script called ``mpif90`` which adds the flags and libraries appropriate
for each type of fortran compiler. So compiling an MPI program usually means simply changing the fortran compiler name
to the MPI script name.

These MPI scripts are built during the MPI install process and are specific to a particular compiler; if your system
has multiple fortran compilers installed then either there will be multiple MPI scripts built, one for each compiler
type, or there will be an environment variable or flag to the MPI script to select which compiler to invoke. See your
system documentation or find an example of a successful MPI program compile command and copy it.

DART use of MPI
~~~~~~~~~~~~~~~

Several DART executables can make use of MPI. To build with MPI, make sure the ``$DART/mkmf/mkmf.template`` has
the correct command for the MPI library you are using. Typically this is mpif90,

.. code::
    
    MPIFC = mpif90
    MPILD = mpif90

but may be ``mpiifort`` if you are using the Intel MPI library.  ``./quickbuild.sh`` will then build DART using MPI.

If you want to build DART without using MPI, run ``./quickbuild.sh nompi``

MPI programs generally need to be started with a shell script called 'mpirun' or 'mpiexec', but they also interact
with any batch control system that might be installed on the cluster or parallel system. Parallel systems with
multiple users generally run some sort of batch system (e.g. LSF, PBS, POE, LoadLeveler, etc). You submit a job
request to this system and it schedules which nodes are assigned to which jobs. Unfortunately the details of this vary
widely from system to system; consult your local web pages or knowledgeable system admin for help here. Generally the
run scripts supplied with DART have generic sections to deal with LSF, PBS, no batch system at all, and sequential
execution, but the details (e.g. the specific queue names, accounting charge codes) will almost certainly have to be
adjusted.

##############
Compiling DART
##############


To build DART executables you will need to 

#. Create an mkmf.template with appropriate compiler and library flags.
#. Choose which model you want to use with DART, and cd into that work directory. 
   For example, if you want to build DART for the lorenz_63 model, ``cd DART/models/lorenz_63/work``
#. Build the DART executables with ``./quickbuild.sh``.

   ``./quickbuild.sh help`` will show the usage options.

   For models:

   .. code-block:: text
   
    ./quickbuild.sh help

    Usage:
     quickbuild.sh                     : build everything
     quickbuild.sh clean               : clean the build
     quickbuild.sh help                : print help message
   
     quickbuild.sh [mpi/nompi/mpif08] [program]   : optional arguments
                                                    [mpi]     build with mpi (default)
                                                    [nompi]   build without mpi
                                                    [mpif08]  build with mpi using mpi_f08
                                                    [program] build a single program
  
     Example 1. Build filter without mpi:
              quickbuild.sh nompi filter
   
     Example 2. Build perfect_model_obs with mpi
              quickbuild.sh perfect_model_obs
   
     Example 3. Build perfect_model_obs with mpi using the mpi_f08 bindings
              quickbuild.sh mpif08 perfect_model_obs


   For observation converters:


   .. code-block:: text

    ./quickbuild.sh help

    Usage:   
     quickbuild.sh               : build everything
     quickbuild.sh clean         : clean the build
     quickbuild.sh help          : print help message
      
     quickbuild.sh [program]     : build a single program



For more detail on how quickbuild.sh works see :ref:`DART build system`.

Customizing the mkmf.template file
=================================================

mkmf.templates for different compilers/architectures can be found
in the ``DART/build_templates`` directory and have names with
extensions that identify the compiler, the architecture, or both. This
is how you inform the build process of the specifics of your system.
**Our intent is that you copy one that is similar to your system into** 
``DART/build_templates/mkmf.template`` **and customize it.**


For the discussion that follows, knowledge of the contents of one of these
templates (e.g. ``DART/build_templates/mkmf.template.intel.linux``)
is needed. Note that only the relevant lines of the file are shown here. The
first portion of the file is a large comment block that provides
valuable advice on how to customize the *mkmf* template file if needed.

.. code-block:: text

   MPIFC = mpif90
   MPILD = mpif90
   FC = ifort
   LD = ifort
   NETCDF = /usr/local
   INCS = -I$(NETCDF)/include
   LIBS = -L$(NETCDF)/lib -lnetcdf -lnetcdff
   FFLAGS = -O2 $(INCS)
   LDFLAGS = $(FFLAGS) $(LIBS)


+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| FC      | the Fortran compiler                                                                                                                                                                                                             |
+=========+==================================================================================================================================================================================================================================+
| LD      | the name of the loader; typically, the same as the Fortran compiler                                                                                                                                                              |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| MPIFC   | the MPI Fortran compiler; see the :doc:`DART MPI introduction <mpi_intro>` for more info                                                                                                                                         |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| MPILD   | the MPI loader; see the :doc:`DART MPI introduction <mpi_intro>` for more info                                                                                                                                                   |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| NETCDF  | the location of your root netCDF installation, which is assumed to contain netcdf.mod and typesizes.mod in the include subdirectory. Note that the value of the NETCDF variable will be used by the “INCS” and “LIBS” variables. |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| INCS    | the includes passed to the compiler during compilation. Note you may need to change this if your netCDF includes netcdf.mod and typesizes.mod are not in the standard location under the include subdirectory of NETCDF.         |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| LIBS    | the libraries passed to “FC” (or “MPIFC”) during compilation. Note you may need to change this if the netCDF libraries libnetcdf and libnetcdff are not in the standard location under the “lib” subdirectory of NETCDF.         |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| FFLAGS  | the Fortran flags passed to “FC” (or “MPIFC”) during compilation. There are often flags used for optimized code versus debugging code. See your particular compiler’s documentation for more information.                        |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| LDFLAGS | the linker flags passed to LD during compilation. See your particular linker’s documentation for more information.                                                                                                               |
+---------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Building the Lorenz_63 DART project
===================================

In order to get started with DART, here we use the Lorenz 63 model,
which is a simple ODE model with only three variables. DART supports
models with many orders of magnitude more variables than three, but if
you can compile and run the DART code for any ONE of the models, you
should be able to compile and run DART for ANY of the models. For
time-dependent filtering known as **cycling**, where observations are
iteratively assimilated at multiple time steps, DART requires the
ability to move the model state forward in time. For low-order models,
this may be possible with a Fortran function call, but for higher-order
models, this is typically done outside of DART’s execution control.
However, the assimilation itself is conducted the same way for **all**
models. For this reason, here we focus solely on the Lorenz 63 model. If
so desired, see :doc:`The Lorenz 63 model: what is it and why should we
care? <lorenz-63-model>` for more information on this simple yet
surprisingly relevant model. See :doc:`A high-level workflow of DA in
DART <high-level-da-workflows>` for further information regarding the DART
workflow if you prefer to do so before building the code.

There are seven separate, stand-alone programs that are typically
necessary for the end-to-end execution of a DART experiment; see below
or the :doc:`What is DART? <what-is-dart>` section for more information on
these programs and their interactions. 

.. note:: some higher-order models have many more than seven programs; for
          example, the Weather Research and Forecasting (WRF) model,
          which is run operationally around the world to predict regional
          weather, has 28 separate programs. Nonetheless, each of these
          programs are built the same way.

To build all DART programs for Lorenz_63:

.. code-block:: text

    cd DART/models/lorenz_63/work
    ./quickbuild.sh


The DART executables are built in the ``work`` directory.
If the build is successful, you will see the following seven programs
in your lorenz_63 work directory.


+--------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Program                                                                                                                  | Purpose                                                                                                                                                                                                                                                                                                         |
+==========================================================================================================================+=================================================================================================================================================================================================================================================================================================================+
|`preprocess   <../assimilation_code/programs/preprocess/preprocess.html>`__                                               | creates custom source code for just the observations of interest                                                                                                                                                                                                                                                |
+--------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|`create_obs_sequence <../assimilation_code/programs/create_obs_sequence/create_obs_sequence.html>`__                      | specify a (set) of observation characteristics taken by a particular (set of) instruments                                                                                                                                                                                                                       |
+--------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|`create_fixed_network_seq <../assimilation_code/programs/create_fixed_network_seq/create_fixed_network_seq.html>`__       | specify the temporal attributes of the observation sets                                                                                                                                                                                                                                                         |
+--------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|`perfect_model_obs <../assimilation_code/programs/perfect_model_obs/perfect_model_obs.html>`__                            | spinup and generate “true state” for synthetic observation experiments                                                                                                                                                                                                                                          |
+--------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|`filter <../assimilation_code/programs/filter/filter.html>`__                                                             | perform data assimilation analysis                                                                                                                                                                                                                                                                              |
+--------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|`obs_diag <../assimilation_code/programs/obs_diag/threed_sphere/obs_diag.html>`__                                         | creates observation-space diagnostic files in netCDF format to support visualization and quantification.                                                                                                                                                                                                        |
+--------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|`obs_sequence_tool <../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html>`__                            | manipulates observation sequence files. This tool is not generally required (particularly for low-order models) but can be used to combine observation sequences or convert from ASCII to binary or vice-versa. Since this is a rather specialized routine, we will not cover its use further in this document. |
+--------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


.. note:: The most common problem is that the netCDF libraries and/or include
          files were not found in the specified location(s). The second most
          common problem is that the netCDF libraries were built with a
          different compiler than the one used for DART. Find (or compile) a 
          compatible netCDF library, edit the ``DART/build_templates/mkmf.template``
          to point to the correct locations of the includes and library files,
          then run ``./quickbuild.sh`` again.

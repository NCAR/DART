##############
Compiling DART
##############

Now that the DART code has been downloaded and the prerequisites have
been verified, you can now begin building and verifying the DART
installation.

Customizing the build scripts — overview
========================================

DART executable programs are constructed using two tools: *mkmf*, and
*make*. The *make* utility is a very commonly used tool that requires a
user-defined input file (a ``Makefile``) that records dependencies
between different source files. *make* then performs actions to the
source hierarchy, in order of dependence, when one or more of the source
files is modified. *mkmf* is a *perl* script that generates a *make*
input file (named *Makefile*) and an example namelist
``input.nml.<program>_default`` with default values.

*mkmf* (think *"make makefile"*) requires two separate input files. The
first is a template file which specifies the commands required for a
specific Fortran90 compiler and may also contain pointers to directories
containing pre- compiled utilities required by the DART system. **This
template file will need to be modified to reflect your system as
detailed in the next section**.

The second input file is a ``path_names`` file which is supplied by DART
and can be used without modification. An *mkmf* command is executed
which uses the ``path_names`` file and the mkmf template file to produce
a ``Makefile`` which is subsequently used by the standard *make*
utility.

Shell scripts that execute the *mkmf* command for all standard DART
executables are provided with the standard DART distribution. For more
information on the `mkmf <https://github.com/NOAA-GFDL/mkmf>`__ tool
please see the `mkmf
documentation <https://extranet.gfdl.noaa.gov/~vb/mkmf.html>`__.

Building and Customizing the ‘mkmf.template’ file
=================================================

A series of templates for different compilers/architectures can be found
in the ``DARTHOME/build_templates`` directory and have names with
extensions that identify the compiler, the architecture, or both. This
is how you inform the build process of the specifics of your system.
**Our intent is that you copy one that is similar to your system into** 
``DARTHOME/build_templates/mkmf.template`` **and customize it.**

For the discussion that follows, knowledge of the contents of one of these
templates (e.g. ``DARTHOME/build_templates/mkmf.template.intel.linux``)
is needed. Note that only the LAST lines of the file are shown here. The
first portion of the file is a large comment block that provides
valuable advice on how to customize the *mkmf* template file if needed.

.. code-block:: bash

   MPIFC = mpif90
   MPILD = mpif90
   FC = ifort
   LD = ifort
   NETCDF = /usr/local
   INCS = -I$(NETCDF)/include
   LIBS = -L$(NETCDF)/lib -lnetcdf -lnetcdff
   FFLAGS = -O2 $(INCS)
   LDFLAGS = $(FFLAGS) $(LIBS)

+-------------------------------------------+--------------------------+
| variable                                  | value                    |
+===========================================+==========================+
| FC                                        | the Fortran compiler     |
+-------------------------------------------+--------------------------+
| LD                                        | the name of the loader;  |
|                                           | typically, the same as   |
|                                           | the Fortran compiler     |
+-------------------------------------------+--------------------------+
| MPIFC                                     | the MPI Fortran          |
|                                           | compiler; see `the DART  |
|                                           | MPI                      |
|                                           | intro                    |
|                                           | duction <dart_mpi.md>`__ |
|                                           | for more info            |
+-------------------------------------------+--------------------------+
| MPILD                                     | the MPI loader; see `the |
|                                           | DART MPI                 |
|                                           | intro                    |
|                                           | duction <dart_mpi.md>`__ |
|                                           | for more info            |
+-------------------------------------------+--------------------------+
| NETCDF                                    | the location of your     |
|                                           | root netCDF              |
|                                           | installation, which is   |
|                                           | assumed to contain       |
|                                           | ``netcdf.mod`` and       |
|                                           | ``typesizes.mod`` in the |
|                                           | include subdirectory.    |
|                                           | Note that the value of   |
|                                           | the NETCDF variable will |
|                                           | be used by the "INCS"    |
|                                           | and "LIBS" variables.    |
+-------------------------------------------+--------------------------+
| INCS                                      | the includes passed to   |
|                                           | the compiler during      |
|                                           | compilation. Note you    |
|                                           | may need to change this  |
|                                           | if your netCDF includes  |
|                                           | ``netcdf.mod`` and       |
|                                           | ``typesizes.mod`` are    |
|                                           | not in the standard      |
|                                           | location under the       |
|                                           | ``include`` subdirectory |
|                                           | of NETCDF.               |
+-------------------------------------------+--------------------------+
| LIBS                                      | the libraries passed to  |
|                                           | "FC" (or "MPIFC") during |
|                                           | compilation. Note you    |
|                                           | may need to change this  |
|                                           | if the netCDF libraries  |
|                                           | ``libnetcdf`` and        |
|                                           | ``libnetcdff`` are not   |
|                                           | in the standard location |
|                                           | under the "lib"          |
|                                           | subdirectory of NETCDF.  |
+-------------------------------------------+--------------------------+
| FFLAGS                                    | the Fortran flags passed |
|                                           | to "FC" (or "MPIFC")     |
|                                           | during compilation.      |
|                                           | There are often flags    |
|                                           | used for optimized code  |
|                                           | versus debugging code.   |
|                                           | See your particular      |
|                                           | compiler’s documentation |
|                                           | for more information.    |
+-------------------------------------------+--------------------------+
| LDFLAGS                                   | the linker flags passed  |
|                                           | to *LD* during           |
|                                           | compilation. See your    |
|                                           | particular linker’s      |
|                                           | documentation for more   |
|                                           | information.             |
+-------------------------------------------+--------------------------+

Customizing the path names files
================================

Several ``path_names_*`` files are provided in the "work" directory for
each specific model. In this case, the directory of interest is
``DARTHOME/models/lorenz_63/work`` (see the next section). Since each
model comes with its own set of files, the ``path_names_*`` files
typically need no customization. However, modifying these files will be
required if you wish to add your model to DART. See `How do I run DART
with my model? <#RunWithMyModel>`__ for more information.

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
so desired, see `The Lorenz 63 model: what is it and why should we
care? <#Lorenz63>`__ for more information on this simple yet
surprisingly relevant model. See `A high-level workflow of DA in
DART <#dartWorkflow>`__ for further information regarding the DART
workflow if you prefer to do so before building the code.

There are seven separate, stand-alone programs that are typically
necessary for the end-to-end execution of a DART experiment; see below
or the `What is DART? <#WhatIsDART>`__ section for more information on
these programs and their interactions. All DART programs are compiled
the same way, and each model directory has a directory called ``work``
that has the components necessary to build the executables.

.. note:: some higher-order models have many more than seven programs; for
          example, the Weather Research and Forecasting (WRF) model,
          which is run operationally around the world to predict regional
          weather, has 28 separate programs. Nonetheless, each of these
          programs are built the same way.

The ``quickbuild.csh`` in each directory builds all seven programs
necessary for Lorenz 63. Describing what the ``quickbuild.csh`` script
does is useful for understanding how to get started with DART.

The following shell commands show how to build two of these seven
programs for the lorenz_63 model: *preprocess* and *obs_diag*.
*preprocess* is a special program that needs to be built and run to
automatically generate Fortran code that is used by DART to support a
subset of observations - which are (potentially) different for every
model. Once *preprocess* has been run and the required Fortran code has
been generated, any of the other DART programs may be built in the same
way as *obs_diag* in this example. Thus, the following runs *mkmf* to
make a ``Makefile`` for *preprocess*, makes the *preprocess* program,
runs *preprocess* to generate the Fortran observation code, runs *mkmf*
to make a ``Makefile`` for *obs_diag*, then makes the *obs_diag*
program:

.. code-block:: bash

   $ cd DARTHOME/models/lorenz_63/work
   $ ./mkmf_preprocess
   $ make
   $ ./preprocess
   $ ./mkmf_obs_diag
   $ make

The remaining executables are built in the same fashion as *obs_diag*:
run the particular *mkmf* script to generate a Makefile, then execute
*make* to build the corresponding program.

Currently, DART executables are built in a ``work`` subdirectory under
the directory containing code for the given model. The Lorenz_63 model
has seven ``mkmf_xxxxxx`` files for the following programs:

+-----------------------------------+-----------------------------------+
| Program                           | Purpose                           |
+===================================+===================================+
| `preproces                        | creates custom source code for    |
| s <../../assimilation_code/progra | just the observations of interest |
| ms/preprocess/preprocess.html>`__ |                                   |
+-----------------------------------+-----------------------------------+
| `cre                              | specify a (set) of observation    |
| ate_obs_sequence <../../assimilat | characteristics taken by a        |
| ion_code/programs/create_obs_sequ | particular (set of) instruments   |
| ence/create_obs_sequence.html>`__ |                                   |
+-----------------------------------+-----------------------------------+
| `create_fixed_netwo               | specify the temporal attributes   |
| rk_seq <../../assimilation_code/p | of the observation sets           |
| rograms/create_fixed_network_seq/ |                                   |
| create_fixed_network_seq.html>`__ |                                   |
+-----------------------------------+-----------------------------------+
| `perfect_model_obs <../../assim   | spinup and generate "true state"  |
| ilation_code/programs/perfect_mod | for synthetic observation         |
| el_obs/perfect_model_obs.html>`__ | experiments                       |
+-----------------------------------+-----------------------------------+
| `filter <../../assimilation_cod   | perform data assimilation         |
| e/programs/filter/filter.html>`__ | analysis                          |
+-----------------------------------+-----------------------------------+
| `obs_diag <../../a                | creates observation-space         |
| ssimilation_code/programs/obs_dia | diagnostic files in netCDF format |
| g/threed_sphere/obs_diag.html>`__ | to support visualization and      |
|                                   | quantification.                   |
+-----------------------------------+-----------------------------------+
| `obs_sequence_tool <../../assim   | manipulates observation sequence  |
| ilation_code/programs/obs_sequenc | files. This tool is not generally |
| e_tool/obs_sequence_tool.html>`__ | required (particularly for        |
|                                   | low-order models) but can be used |
|                                   | to combine observation sequences  |
|                                   | or convert from ASCII to binary   |
|                                   | or vice-versa. Since this is a    |
|                                   | rather specialized routine, we    |
|                                   | will not cover its use further in |
|                                   | this document.                    |
+-----------------------------------+-----------------------------------+

As mentioned above, ``quickbuild.csh`` is a script that will build every
executable in the directory. There is an optional argument that will
additionally build the MPI-enabled versions which will not be covered in
this set of instructions. See The DART MPI introduction page for more
information on using DART with MPI.

Running ``quickbuild.csh`` will compile all the executables mentioned
above for the lorenz_63 model:

.. code-block:: bash

   $ cd DARTHOME/models/lorenz_63/work
   $ ./quickbuild.csh

The result (hopefully) is that seven executables now reside in your work
directory.

.. note:: The most common problem is that the netCDF libraries and/or include
          files were not found in the specified location(s). The second most
          common problem is that the netCDF libraries were built with a
          different compiler than the one used for DART. Find (or compile) a 
          compatible netCDF library, edit the ``DARTHOME/build_templates/mkmf.template``
          to point to the correct locations of the includes and library files,
          recreate the ``Makefile``\ s, and try again.
Kodiak
======

DART Kodiak release documentation
---------------------------------

.. attention::

   **Kodiak** was released in June of 2011. Its source code is available via the `DART repository on
   Github <https://github.com/NCAR/DART/tree/Kodiak>`__. This documentation is preserved merely for reference. See the
   `DART homepage <https://dart.ucar.edu/>`__ to learn about the latest release.

DART overview
-------------

The Data Assimilation Research Testbed (DART) is designed to facilitate the combination of assimilation algorithms,
models, and real (or synthetic) observations to allow increased understanding of all three. The DART programs are highly
portable, having been compiled with many Fortran 90 compilers and run on linux compute-servers, linux clusters, OSX
laptops/desktops, SGI Altix clusters, supercomputers running AIX, and more. Read the Customizations section for help in
building on new platforms.

DART employs a modular programming approach to apply an Ensemble Kalman Filter which nudges models toward a state that
is more consistent with information from a set of observations. Models may be swapped in and out, as can different
algorithms in the Ensemble Kalman Filter. The method requires running multiple instances of a model to generate an
ensemble of states. A forward operator appropriate for the type of observation being assimilated is applied to each of
the states to generate the model's estimate of the observation. Comparing these estimates and their uncertainty to the
observation and its uncertainty ultimately results in the adjustments to the model states. There's much more to it,
described in detail in the tutorial directory of the package.

DART diagnostic output includes two netCDF files containing the model states just before the adjustment
(``Prior_Diag.nc``) and just after the adjustment (``Posterior_Diag.nc``) as well as a file ``obs_seq.final`` with the
model estimates of the observations. There is a suite of Matlab® functions that facilitate exploration of the results,
but the netCDF files are inherently portable and contain all the necessary metadata to interpret the contents.

In this document links are available which point to Web-based documentation files and also to the same information in
html files distributed with DART. If you have used subversion to check out a local copy of the DART files you can open
this file in a browser by loading ``DART/doc/html/Kodiak_release.html`` and then use the ``local file`` links to see
other documentation pages without requiring a connection to the internet. If you are looking at this documentation from
the ``www.image.ucar.edu`` web server or you are connected to the internet you can use the ``Website`` links to view
other documentation pages.

Getting started
---------------

What's required
~~~~~~~~~~~~~~~

#. a Fortran 90 compiler
#. a netCDF library including the F90 interfaces
#. the C shell
#. (optional, to run in parallel) an MPI library

DART has been tested on many Fortran compilers and platforms. We don't have any platform-dependent code sections and we
use only the parts of the language that are portable across all the compilers we have access to. We explicitly set the
Fortran 'kind' for all real values and do not rely on autopromotion or other compile-time flags to set the default byte
size for numbers. It is possible that some model-specific interface code from outside sources may have specific compiler
flag requirements; see the documentation for each model. The low-order models and all common portions of the DART code
compile cleanly.

| DART uses the `netCDF <http://www.unidata.ucar.edu/packages/netcdf/>`__ self-describing data format with a particular
  metadata convention to describe output that is used to analyze the results of assimilation experiments. These files
  have the extension ``.nc`` and can be read by a number of standard data analysis tools.
| Since most of the models being used with DART are written in Fortran and run on various UNIX or \*nix platforms, the
  development environment for DART is highly skewed to these machines. We do most of our development on a small linux
  workstation and a mac laptop running OSX 10.x, and we have an extensive test network. (I've never built nor run DART
  on a Windows machine - so I don't even know if it's possible. If you have run it (under Cygwin?) please let me know
  how it went -- I'm curious. Tim - thoar 'at' ucar 'dot ' edu)

What's nice to have
~~~~~~~~~~~~~~~~~~~

-  **ncview**: DART users have used `ncview <http://meteora.ucsd.edu/~pierce/ncview_home_page.html>`__ to create
   graphical displays of output data fields. The 2D rendering is good for 'quick-look' type uses, but I wouldn't want to
   publish with it.
-  **NCO**: The `NCO <http://nco.sourceforge.net>`__ tools are able to perform operations on netCDF files like
   concatenating, slicing, and dicing.
-  **Matlab**\ ®: A set of `Matlab® <http://www.mathworks.com/>`__ scripts designed to produce graphical diagnostics
   from DART netCDF output files are also part of the DART project.
-  **MPI**: The DART system includes an MPI option. MPI stands for 'Message Passing Interface', and is both a library
   and run-time system that enables multiple copies of a single program to run in parallel, exchange data, and combine
   to solve a problem more quickly. DART does **NOT** require MPI to run; the default build scripts do not need nor use
   MPI in any way. However, for larger models with large state vectors and large numbers of observations, the data
   assimilation step will run much faster in parallel, which requires MPI to be installed and used. However, if multiple
   ensembles of your model fit comfortably (in time and memory space) on a single processor, you need read no further
   about MPI.

Types of input
~~~~~~~~~~~~~~

DART programs can require three different types of input. First, some of the DART programs, like those for creating
synthetic observational datasets, require interactive input from the keyboard. For simple cases this interactive input
can be made directly from the keyboard. In more complicated cases a file containing the appropriate keyboard input can
be created and this file can be directed to the standard input of the DART program. Second, many DART programs expect
one or more input files in DART specific formats to be available. For instance, ``perfect_model_obs``, which creates a
synthetic observation set given a particular model and a description of a sequence of observations, requires an input
file that describes this observation sequence. At present, the observation files for DART are in a custom format in
either human-readable ascii or more compact machine-specific binary. Third, many DART modules (including main programs)
make use of the Fortran90 namelist facility to obtain values of certain parameters at run-time. All programs look for a
namelist input file called ``input.nml`` in the directory in which the program is executed. The ``input.nml`` file can
contain a sequence of individual Fortran90 namelists which specify values of particular parameters for modules that
compose the executable program.

Installation
------------

This document outlines the installation of the DART software and the system requirements. The entire installation
process is summarized in the following steps:

#. Determine which F90 compiler is available.
#. Determine the location of the ``netCDF`` library.
#. Download the DART software into the expected source tree.
#. Modify certain DART files to reflect the available F90 compiler and location of the appropriate libraries.
#. Build the executables.

We have tried to make the code as portable as possible, but we do not have access to all compilers on all platforms, so
there are no guarantees. We are interested in your experience building the system, so please email me (Tim Hoar) thoar
'at' ucar 'dot' edu (trying to cut down on the spam).

After the installation, you might want to peruse the following.

-  Running the Lorenz_63 Model.
-  Using the Matlab® diagnostic scripts.
-  A short discussion on bias, filter divergence and covariance inflation.
-  And another one on synthetic observations.

You should *absolutely* run the DART_LAB interactive tutorial (if you have Matlab available) and look at the DART_LAB
presentation slides `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/DART_LAB/DART_LAB.html>`__ or
:doc:`../DART_LAB/DART_LAB` in the ``DART_LAB`` directory, and then take the tutorial in the ``DART/tutorial``
directory.

Requirements: an F90 compiler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DART software has been successfully built on several Linux/x86 platforms with several versions of the `Intel Fortran
Compiler for Linux <http://www.intel.com/software/products/compilers/flin>`__, which (at one point) is/was free for
individual scientific use. Also Intel Fortran for Mac OSX. It has also been built and successfully run with several
versions of each of the following: `Portland Group Fortran Compiler <http://www.pgroup.com>`__, `Lahey Fortran
Compiler <http://www.lahey.com>`__, `Pathscale Fortran Compiler <http://www.pathscale.com>`__, `GNU Fortran 95 Compiler
("gfortran") <http://gcc.gnu.org/fortran>`__, `Absoft Fortran 90/95 Compiler (Mac OSX) <http://www.absoft.com>`__. Since
recompiling the code is a necessity to experiment with different models, there are no binaries to distribute.

DART uses the `netCDF <http://www.unidata.ucar.edu/packages/netcdf/>`__ self-describing data format for the results of
assimilation experiments. These files have the extension ``.nc`` and can be read by a number of standard data analysis
tools. In particular, DART also makes use of the F90 interface to the library which is available through the
``netcdf.mod`` and ``typesizes.mod`` modules. *IMPORTANT*: different compilers create these modules with different
"case" filenames, and sometimes they are not **both** installed into the expected directory. It is required that both
modules be present. The normal place would be in the ``netcdf/include`` directory, as opposed to the ``netcdf/lib``
directory.

If the netCDF library does not exist on your system, you must build it (as well as the F90 interface modules). The
library and instructions for building the library or installing from an RPM may be found at the netCDF home page:
http://www.unidata.ucar.edu/packages/netcdf/ Pay particular attention to the compiler-specific patches that must be
applied for the Intel Fortran Compiler. (Or the PG compiler, for that matter.)

The location of the netCDF library, ``libnetcdf.a``, and the locations of both ``netcdf.mod`` and ``typesizes.mod`` will
be needed by the makefile template, as described in the compiling section. Depending on the netCDF build options, the
Fortran 90 interfaces may be built in a separate library named ``netcdff.a`` and you may need to add ``-lnetcdff`` to
the library flags.

Unpacking the distribution
--------------------------

This release of the `DART source code can be downloaded <https://github.com/NCAR/DART/releases/tag/v7.0.0>`__ as a
compressed zip or tar.gz file. When extracted, the source tree will begin with a directory named ``DART`` and will be
approximately 206.5 Mb. Compiling the code in this tree (as is usually the case) will necessitate much more space.

::


   $ gunzip DART-7.0.0.tar.gz
   $ tar -xvf DART-7.0.0.tar

You should wind up with a directory named ``DART``.

The code tree is very "bushy"; there are many directories of support routines, etc. but only a few directories involved
with the customization and installation of the DART software. If you can compile and run ONE of the low-order models,
you should be able to compile and run ANY of the low-order models. For this reason, we can focus on the Lorenz \`63
model. Subsequently, the only directories with files to be modified to check the installation are: ``DART/mkmf``,
``DART/models/lorenz_63/work``, and ``DART/matlab`` (but only for analysis).

Customizing the build scripts -- overview
-----------------------------------------

DART executable programs are constructed using two tools: ``make`` and ``mkmf``. The ``make`` utility is a very common
piece of software that requires a user-defined input file that records dependencies between different source files.
``make`` then performs a hierarchy of actions when one or more of the source files is modified. The ``mkmf`` utility is
a custom preprocessor that generates a ``make`` input file (named ``Makefile``) and an example namelist
*input.nml.\ program\ \_default* with the default values. The ``Makefile`` is designed specifically to work with
object-oriented Fortran90 (and other languages) for systems like DART.

``mkmf`` requires two separate input files. The first is a \`template' file which specifies details of the commands
required for a specific Fortran90 compiler and may also contain pointers to directories containing pre-compiled
utilities required by the DART system. **This template file will need to be modified to reflect your system**. The
second input file is a \`path_names' file which includes a complete list of the locations (either relative or absolute)
of all Fortran90 source files that are required to produce a particular DART program. Each 'path_names' file must
contain a path for exactly one Fortran90 file containing a main program, but may contain any number of additional paths
pointing to files containing Fortran90 modules. An ``mkmf`` command is executed which uses the 'path_names' file and the
mkmf template file to produce a ``Makefile`` which is subsequently used by the standard ``make`` utility.

Shell scripts that execute the mkmf command for all standard DART executables are provided as part of the standard DART
software. For more information on ``mkmf`` see `the FMS mkmf
description <http://www.gfdl.gov/fms/pubrel/j/atm_dycores/doc/dycore_public_manual.html#mkmf>`__.

One of the benefits of using ``mkmf`` is that it also creates an example namelist file for each program. The example
namelist is called *input.nml.\ program\ \_default*, so as not to clash with any exising ``input.nml`` that may exist in
that directory.

Building and customizing the 'mkmf.template' file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A series of templates for different compilers/architectures exists in the ``DART/mkmf/`` directory and have names with
extensions that identify the compiler, the architecture, or both. This is how you inform the build process of the
specifics of your system. Our intent is that you copy one that is similar to your system into ``mkmf.template`` and
customize it. For the discussion that follows, knowledge of the contents of one of these templates (i.e.
``mkmf.template.gfortran``) is needed. Note that only the LAST lines are shown here, the head of the file is just a big
comment (worth reading, btw).

::


   ...
   MPIFC = mpif90
   MPILD = mpif90
   FC = gfortran
   LD = gfortran
   NETCDF = /usr/local
   INCS = ${NETCDF}/include
   FFLAGS = -O2 -I$(INCS)
   LIBS = -L${NETCDF}/lib -lnetcdf
   LDFLAGS = -I$(INCS) $(LIBS)

Essentially, each of the lines defines some part of the resulting ``Makefile``. Since ``make`` is particularly good at
sorting out dependencies, the order of these lines really doesn't make any difference. The ``FC = gfortran`` line
ultimately defines the Fortran90 compiler to use, etc. The lines which are most likely to need site-specific changes
start with ``FFLAGS`` and ``NETCDF``, which indicate where to look for the netCDF F90 modules and the location of the
netCDF library and modules.

If you have MPI installed on your system ``MPIFC, MPILD`` dictate which compiler will be used in that instance. If you
do not have MPI, these variables are of no consequence.

Netcdf
^^^^^^

| Modifying the ``NETCDF`` value should be relatively straightforward.
| Change the string to reflect the location of your netCDF installation containing ``netcdf.mod`` and ``typesizes.mod``.
  The value of the ``NETCDF`` variable will be used by the ``FFLAGS, LIBS,`` and ``LDFLAGS`` variables.

FFLAGS
^^^^^^

Each compiler has different compile flags, so there is really no way to exhaustively cover this other than to say the
templates as we supply them should work -- depending on the location of your netCDF. The low-order models can be
compiled without a ``-r8`` switch, but the ``bgrid_solo`` model cannot.

Libs
^^^^

The Fortran 90 interfaces may be part of the default ``netcdf.a`` library and ``-lnetcdf`` is all you need. However it
is also common for the Fortran 90 interfaces to be built in a separate library named ``netcdff.a``. In that case you
will need ``-lnetcdf`` and also ``-lnetcdff`` on the **LIBS** line. This is a build-time option when the netCDF
libraries are compiled so it varies from site to site.

Customizing the 'path_names_*' file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several ``path_names_*`` files are provided in the ``work`` directory for each specific model, in this case:
``DART/models/lorenz_63/work``. Since each model comes with its own set of files, the ``path_names_*`` files need no
customization.

Building the Lorenz_63 DART project
-----------------------------------

DART executables are constructed in a ``work`` subdirectory under the directory containing code for the given model.
From the top-level DART directory change to the L63 work directory and list the contents:

::


   $ cd DART/models/lorenz_63/work
   $ ls -1

With the result:

::


   Posterior_Diag.nc
   Prior_Diag.nc
   True_State.nc
   filter_ics
   filter_restart
   input.nml
   mkmf_create_fixed_network_seq
   mkmf_create_obs_sequence
   mkmf_filter
   mkmf_obs_diag
   mkmf_obs_sequence_tool
   mkmf_perfect_model_obs
   mkmf_preprocess
   mkmf_restart_file_tool
   mkmf_wakeup_filter
   obs_seq.final
   obs_seq.in
   obs_seq.out
   obs_seq.out.average
   obs_seq.out.x
   obs_seq.out.xy
   obs_seq.out.xyz
   obs_seq.out.z
   path_names_create_fixed_network_seq
   path_names_create_obs_sequence
   path_names_filter
   path_names_obs_diag
   path_names_obs_sequence_tool
   path_names_perfect_model_obs
   path_names_preprocess
   path_names_restart_file_tool
   path_names_wakeup_filter
   perfect_ics
   perfect_restart
   quickbuild.csh
   set_def.out
   workshop_setup.csh

In all the ``work`` directories there will be a ``quickbuild.csh`` script that builds or rebuilds the executables. The
following instructions do this work by hand to introduce you to the individual steps, but in practice running quickbuild
will be the normal way to do the compiles.

There are nine ``mkmf_``\ *xxxxxx* files for the programs

#. ``preprocess``,
#. ``create_obs_sequence``,
#. ``create_fixed_network_seq``,
#. ``perfect_model_obs``,
#. ``filter``,
#. ``wakeup_filter``,
#. ``obs_sequence_tool``, and
#. ``restart_file_tool``, and
#. ``obs_diag``,

along with the corresponding ``path_names_``\ *xxxxxx* files. There are also files that contain initial conditions,
netCDF output, and several observation sequence files, all of which will be discussed later. You can examine the
contents of one of the ``path_names_``\ *xxxxxx* files, for instance ``path_names_filter``, to see a list of the
relative paths of all files that contain Fortran90 modules required for the program ``filter`` for the L63 model. All of
these paths are relative to your ``DART`` directory. The first path is the main program (``filter.f90``) and is followed
by all the Fortran90 modules used by this program (after preprocessing).

The ``mkmf_``\ *xxxxxx* scripts are cryptic but should not need to be modified -- as long as you do not restructure the
code tree (by moving directories, for example). The function of the ``mkmf_``\ *xxxxxx* script is to generate a
``Makefile`` and an *input.nml.\ program\ \_default* file. It does not do the compile; ``make`` does that:

::


   $ csh mkmf_preprocess
   $ make

The first command generates an appropriate ``Makefile`` and the ``input.nml.preprocess_default`` file. The second
command results in the compilation of a series of Fortran90 modules which ultimately produces an executable file:
``preprocess``. Should you need to make any changes to the ``DART/mkmf/mkmf.template``, you will need to regenerate the
``Makefile``.

The ``preprocess`` program actually builds source code to be used by all the remaining modules. It is **imperative** to
actually **run** ``preprocess`` before building the remaining executables. This is how the same code can assimilate
state vector 'observations' for the Lorenz_63 model and real radar reflectivities for WRF without needing to specify a
set of radar operators for the Lorenz_63 model!

``preprocess`` reads the ``&preprocess_nml`` namelist to determine what observations and operators to incorporate. For
this exercise, we will use the values in ``input.nml``. ``preprocess`` is designed to abort if the files it is supposed
to build already exist. For this reason, it is necessary to remove a couple files (if they exist) before you run the
preprocessor. (The ``quickbuild.csh`` script will do this for you automatically.)

.. container:: unix

   ::

      $ \rm -f ../../obs_def/obs_def_mod.f90
      $ \rm -f ../../obs_kind/obs_kind_mod.f90
      $ ./preprocess
      $ ls -l ../../obs_def/obs_def_mod.f90
      $ ls -l ../../obs_kind/obs_kind_mod.f90

This created ``../../obs_def/obs_def_mod.f90`` from ``../../obs_kind/DEFAULT_obs_kind_mod.F90`` and several other
modules. ``../../obs_kind/obs_kind_mod.f90`` was created similarly. Now we can build the rest of the project.

A series of object files for each module compiled will also be left in the work directory, as some of these are
undoubtedly needed by the build of the other DART components. You can proceed to create the other programs needed to
work with L63 in DART as follows:

::


   $ csh mkmf_create_obs_sequence
   $ make
   $ csh mkmf_create_fixed_network_seq
   $ make
   $ csh mkmf_perfect_model_obs
   $ make
   $ csh mkmf_filter
   $ make
   $ csh mkmf_obs_diag
   $ make

The result (hopefully) is that six executables now reside in your work directory. The most common problem is that the
netCDF libraries and include files (particularly ``typesizes.mod``) are not found. Edit the ``DART/mkmf/mkmf.template``,
recreate the ``Makefile``, and try again.

+------------------------------+--------------------------------------------------------------------------------------+
| program                      | purpose                                                                              |
+==============================+======================================================================================+
| ``preprocess``               | creates custom source code for just the observation types of interest                |
+------------------------------+--------------------------------------------------------------------------------------+
| ``create_obs_sequence``      | specify a (set) of observation characteristics taken by a particular (set of)        |
|                              | instruments                                                                          |
+------------------------------+--------------------------------------------------------------------------------------+
| ``create_fixed_network_seq`` | repeat a set of observations through time to simulate observing networks where       |
|                              | observations are taken in the same location at regular (or irregular) intervals      |
+------------------------------+--------------------------------------------------------------------------------------+
| ``perfect_model_obs``        | generate "true state" for synthetic observation experiments. Can also be used to     |
|                              | 'spin up' a model by running it for a long time.                                     |
+------------------------------+--------------------------------------------------------------------------------------+
| ``filter``                   | does the assimilation                                                                |
+------------------------------+--------------------------------------------------------------------------------------+
| ``obs_diag``                 | creates observation-space diagnostic files to be explored by the Matlab® scripts.    |
+------------------------------+--------------------------------------------------------------------------------------+
| ``obs_sequence_tool``        | manipulates observation sequence files. It is not generally needed (particularly for |
|                              | low-order models) but can be used to combine observation sequences or convert from   |
|                              | ASCII to binary or vice-versa. We will not cover its use in this document.           |
+------------------------------+--------------------------------------------------------------------------------------+
| ``restart_file_tool``        | manipulates the initial condition and restart files. We're going to ignore this one  |
|                              | here.                                                                                |
+------------------------------+--------------------------------------------------------------------------------------+
| ``wakeup_filter``            | is only needed for MPI applications. We're starting at the beginning here, so we're  |
|                              | going to ignore this one, too.                                                       |
+------------------------------+--------------------------------------------------------------------------------------+

Running Lorenz_63
-----------------

This initial sequence of exercises includes detailed instructions on how to work with the DART code and allows
investigation of the basic features of one of the most famous dynamical systems, the 3-variable Lorenz-63 model. The
remarkable complexity of this simple model will also be used as a case study to introduce a number of features of a
simple ensemble filter data assimilation system. To perform a synthetic observation assimilation experiment for the L63
model, the following steps must be performed (an overview of the process is given first, followed by detailed procedures
for each step):

Experiment overview
-------------------

#. Integrate the L63 model for a long time
   starting from arbitrary initial conditions to generate a model state that lies on the attractor. The ergodic nature
   of the L63 system means a 'lengthy' integration always converges to some point on the computer's finite precision
   representation of the model's attractor.
#. Generate a set of ensemble initial conditions
   from which to start an assimilation. Since L63 is ergodic, the ensemble members can be designed to look like random
   samples from the model's 'climatological distribution'. To generate an ensemble member, very small perturbations can
   be introduced to the state on the attractor generated by step 1. This perturbed state can then be integrated for a
   very long time until all memory of its initial condition can be viewed as forgotten. Any number of ensemble initial
   conditions can be generated by repeating this procedure.
#. Simulate a particular observing system
   by first creating an 'observation set definition' and then creating an 'observation sequence'. The 'observation set
   definition' describes the instrumental characteristics of the observations and the 'observation sequence' defines the
   temporal sequence of the observations.
#. Populate the 'observation sequence' with 'perfect' observations
   by integrating the model and using the information in the 'observation sequence' file to create simulated
   observations. This entails operating on the model state at the time of the observation with an appropriate forward
   operator (a function that operates on the model state vector to produce the expected value of the particular
   observation) and then adding a random sample from the observation error distribution specified in the observation set
   definition. At the same time, diagnostic output about the 'true' state trajectory can be created.
#. Assimilate the synthetic observations
   by running the filter; diagnostic output is generated.

1. Integrate the L63 model for a 'long' time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``perfect_model_obs`` integrates the model for all the times specified in the 'observation sequence definition' file. To
this end, begin by creating an 'observation sequence definition' file that spans a long time. Creating an 'observation
sequence definition' file is a two-step procedure involving ``create_obs_sequence`` followed by
``create_fixed_network_seq``. After they are both run, it is necessary to integrate the model with
``perfect_model_obs``.

1.1 Create an observation set definition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``create_obs_sequence`` creates an observation set definition, the time-independent part of an observation sequence. An
observation set definition file only contains the ``location, type,`` and ``observational error characteristics``
(normally just the diagonal observational error variance) for a related set of observations. There are no actual
observations, nor are there any times associated with the definition. For spin-up, we are only interested in integrating
the L63 model, not in generating any particular synthetic observations. Begin by creating a minimal observation set
definition.

In general, for the low-order models, only a single observation set need be defined. Next, the number of individual
scalar observations (like a single surface pressure observation) in the set is needed. To spin-up an initial condition
for the L63 model, only a single observation is needed. Next, the error variance for this observation must be entered.
Since we do not need (nor want) this observation to have any impact on an assimilation (it will only be used for
spinning up the model and the ensemble), enter a very large value for the error variance. An observation with a very
large error variance has essentially no impact on deterministic filter assimilations like the default variety
implemented in DART. Finally, the location and type of the observation need to be defined. For all types of models, the
most elementary form of synthetic observations are called 'identity' observations. These observations are generated
simply by adding a random sample from a specified observational error distribution directly to the value of one of the
state variables. This defines the observation as being an identity observation of the first state variable in the L63
model. The program will respond by terminating after generating a file (generally named ``set_def.out``) that defines
the single identity observation of the first state variable of the L63 model. The following is a screenshot (much of the
verbose logging has been left off for clarity), the user input looks *like this*.

.. container:: unix

   ::

      [unixprompt]$ ./create_obs_sequence
       Starting program create_obs_sequence
       Initializing the utilities module.
       Trying to log to unit   10
       Trying to open file dart_log.out
       
       Registering module :
       $url: http://squish/DART/trunk/utilities/utilities_mod.f90 $
       $revision: 2713 $
       $date: 2007-03-25 22:09:04 -0600 (Sun, 25 Mar 2007) $
       Registration complete.

       &UTILITIES_NML
       TERMLEVEL= 2,LOGFILENAME=dart_log.out                                          
                                                                                  
       /
       
       Registering module :
       $url: http://squish/DART/trunk/obs_sequence/create_obs_sequence.f90 $
       $revision: 2713 $
       $date: 2007-03-25 22:09:04 -0600 (Sun, 25 Mar 2007) $
       Registration complete.

       { ... }

       Input upper bound on number of observations in sequence
      10
       
       Input number of copies of data (0 for just a definition)
      0

       Input number of quality control values per field (0 or greater)
      0

       input a -1 if there are no more obs 
      0

       Registering module :
       $url: http://squish/DART/trunk/obs_def/DEFAULT_obs_def_mod.F90 $
       $revision: 2820 $
       $date: 2007-04-09 10:37:47 -0600 (Mon, 09 Apr 2007) $
       Registration complete.
       
       
       Registering module :
       $url: http://squish/DART/trunk/obs_kind/DEFAULT_obs_kind_mod.F90 $
       $revision: 2822 $
       $date: 2007-04-09 10:39:08 -0600 (Mon, 09 Apr 2007) $
       Registration complete.
       
       ------------------------------------------------------
       
       initialize_module obs_kind_nml values are
       
       -------------- ASSIMILATE_THESE_OBS_TYPES --------------
       RAW_STATE_VARIABLE
       -------------- EVALUATE_THESE_OBS_TYPES --------------
       ------------------------------------------------------
       
            Input -1 * state variable index for identity observations
            OR input the name of the observation kind from table below:
            OR input the integer index, BUT see documentation...
              1 RAW_STATE_VARIABLE

      -1

       input time in days and seconds
      1 0

       Input error variance for this observation definition
      1000000

       input a -1 if there are no more obs 
      -1

       Input filename for sequence (  set_def.out   usually works well)
       set_def.out 
       write_obs_seq  opening formatted file set_def.out
       write_obs_seq  closed file set_def.out

1.2 Create an observation sequence definition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``create_fixed_network_seq`` creates an 'observation sequence definition' by extending the 'observation set definition'
with the temporal attributes of the observations.

The first input is the name of the file created in the previous step, i.e. the name of the observation set definition
that you've just created. It is possible to create sequences in which the observation sets are observed at regular
intervals or irregularly in time. Here, all we need is a sequence that takes observations over a long period of time -
indicated by entering a 1. Although the L63 system normally is defined as having a non-dimensional time step, the DART
system arbitrarily defines the model timestep as being 3600 seconds. If we declare that we have one observation per day
for 1000 days, we create an observation sequence definition spanning 24000 'model' timesteps; sufficient to spin-up the
model onto the attractor. Finally, enter a name for the 'observation sequence definition' file. Note again: there are no
observation values present in this file. Just an observation type, location, time and the error characteristics. We are
going to populate the observation sequence with the ``perfect_model_obs`` program.

.. container:: unix

   ::

      [unixprompt]$ ./create_fixed_network_seq

       ...

       Registering module :
       $url: http://squish/DART/trunk/obs_sequence/obs_sequence_mod.f90 $
       $revision: 2749 $
       $date: 2007-03-30 15:07:33 -0600 (Fri, 30 Mar 2007) $
       Registration complete.
       
       static_init_obs_sequence obs_sequence_nml values are
       &OBS_SEQUENCE_NML
       WRITE_BINARY_OBS_SEQUENCE =  F,
       /
       Input filename for network definition sequence (usually  set_def.out  )
      set_def.out

       ...

       To input a regularly repeating time sequence enter 1
       To enter an irregular list of times enter 2
      1
       Input number of observations in sequence
      1000
       Input time of initial ob in sequence in days and seconds
      1, 0
       Input period of obs in days and seconds
      1, 0
                 1
                 2
                 3
      ...
               997
               998
               999
              1000
      What is output file name for sequence (  obs_seq.in   is recommended )
      obs_seq.in
       write_obs_seq  opening formatted file obs_seq.in
       write_obs_seq closed file obs_seq.in

1.3 Initialize the model onto the attractor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``perfect_model_obs`` can now advance the arbitrary initial state for 24,000 timesteps to move it onto the attractor.

``perfect_model_obs`` uses the Fortran90 namelist input mechanism instead of (admittedly gory, but temporary)
interactive input. All of the DART software expects the namelists to found in a file called ``input.nml``. When you
built the executable, an example namelist was created ``input.nml.perfect_model_obs_default`` that contains all of the
namelist input for the executable. If you followed the example, each namelist was saved to a unique name. We must now
rename and edit the namelist file for ``perfect_model_obs``. Copy ``input.nml.perfect_model_obs_default`` to
``input.nml`` and edit it to look like the following: (just worry about the highlighted stuff - and whitespace doesn't
matter)

::


   $ cp input.nml.perfect_model_obs_default
   $ input.nml

.. container:: routineIndent1

   ::

      &perfect_model_obs_nml
         start_from_restart    = .false.,
         output_restart        = .true.,
         async                 = 0,
         init_time_days        = 0,
         init_time_seconds     = 0,
         first_obs_days        = -1,
         first_obs_seconds     = -1,
         last_obs_days         = -1,
         last_obs_seconds      = -1,
         output_interval       = 1,
         restart_in_file_name  = "perfect_ics",
         restart_out_file_name = "perfect_restart",
         obs_seq_in_file_name  = "obs_seq.in",
         obs_seq_out_file_name = "obs_seq.out",
         adv_ens_command       = "./advance_ens.csh"  /

      &ensemble_manager_nml
         single_restart_file_in  = .true.,
         single_restart_file_out = .true.,
         perturbation_amplitude  = 0.2  /

      &assim_tools_nml
         filter_kind                     = 1,
         cutoff                          = 0.2,
         sort_obs_inc                    = .false.,
         spread_restoration              = .false.,
         sampling_error_correction       = .false.,
         adaptive_localization_threshold = -1,
         print_every_nth_obs             = 0  /

      &cov_cutoff_nml
         select_localization = 1  /

      &reg_factor_nml
         select_regression    = 1,
         input_reg_file       = "time_mean_reg",
         save_reg_diagnostics = .false.,
         reg_diagnostics_file = "reg_diagnostics"  /

      &obs_sequence_nml
         write_binary_obs_sequence = .false.  /

      &obs_kind_nml
         assimilate_these_obs_types = 'RAW_STATE_VARIABLE'  /

      &assim_model_nml
         write_binary_restart_files = .true. /

      &model_nml
         sigma  = 10.0,
         r      = 28.0,
         b      = 2.6666666666667,
         deltat = 0.01,
         time_step_days = 0,
         time_step_seconds = 3600  /

      &utilities_nml
         TERMLEVEL = 1,
         logfilename = 'dart_log.out'  /

For the moment, only two namelists warrant explanation. Each namelists is covered in detail in the html files
accompanying the source code for the module.

perfect_model_obs_nml
~~~~~~~~~~~~~~~~~~~~~

+---------------------------+-----------------------------------------------------------------------------------------+
| namelist variable         | description                                                                             |
+===========================+=========================================================================================+
| ``start_from_restart``    | When set to 'false', ``perfect_model_obs`` generates an arbitrary initial condition     |
|                           | (which cannot be guaranteed to be on the L63 attractor). When set to 'true', a restart  |
|                           | file (specified by ``restart_in_file_name``) is read.                                   |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``output_restart``        | When set to 'true', ``perfect_model_obs`` will record the model state at the end of     |
|                           | this integration in the file named by ``restart_out_file_name``.                        |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``async``                 | The lorenz_63 model is advanced through a subroutine call - indicated by async = 0.     |
|                           | There is no other valid value for this model.                                           |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``init_time_``\ *xxxx*    | the start time of the integration.                                                      |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``first_obs_``\ *xxxx*    | the time of the first observation of interest. While not needed in this example, you    |
|                           | can skip observations if you want to. A value of -1 indicates to start at the           |
|                           | beginning.                                                                              |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``last_obs_``\ *xxxx*     | the time of the last observation of interest. While not needed in this example, you do  |
|                           | not have to assimilate all the way to the end of the observation sequence file. A value |
|                           | of -1 indicates to use all the observations.                                            |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``output_interval``       | interval at which to save the model state (in True_State.nc).                           |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``restart_in_file_name``  | is ignored when 'start_from_restart' is 'false'.                                        |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``restart_out_file_name`` | if ``output_restart`` is 'true', this specifies the name of the file containing the     |
|                           | model state at the end of the integration.                                              |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``obs_seq_in_file_name``  | specifies the file name that results from running ``create_fixed_network_seq``, i.e.    |
|                           | the 'observation sequence definition' file.                                             |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``obs_seq_out_file_name`` | specifies the output file name containing the 'observation sequence', finally populated |
|                           | with (perfect?) 'observations'.                                                         |
+---------------------------+-----------------------------------------------------------------------------------------+
| ``advance_ens_command``   | specifies the shell commands or script to execute when async /= 0.                      |
+---------------------------+-----------------------------------------------------------------------------------------+

utilities_nml
~~~~~~~~~~~~~

+-------------------+-------------------------------------------------------------------------------------------------+
| namelist variable | description                                                                                     |
+===================+=================================================================================================+
| ``TERMLEVEL``     | When set to '1' the programs terminate when a 'warning' is generated. When set to '2' the       |
|                   | programs terminate only with 'fatal' errors.                                                    |
+-------------------+-------------------------------------------------------------------------------------------------+
| ``logfilename``   | Run-time diagnostics are saved to this file. This namelist is used by all programs, so the file |
|                   | is opened in APPEND mode. Subsequent executions cause this file to grow.                        |
+-------------------+-------------------------------------------------------------------------------------------------+

Executing ``perfect_model_obs`` will integrate the model 24,000 steps and output the resulting state in the file
``perfect_restart``. Interested parties can check the spinup in the ``True_State.nc`` file.

::


   $ ./perfect_model_obs

2. Generate a set of ensemble initial conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The set of initial conditions for a 'perfect model' experiment is created in several steps. 1) Starting from the spun-up
state of the model (available in ``perfect_restart``), run ``perfect_model_obs`` to generate the 'true state' of the
experiment and a corresponding set of observations. 2) Feed the same initial spun-up state and resulting observations
into ``filter``.

The first step is achieved by changing a perfect_model_obs namelist parameter, copying ``perfect_restart`` to
``perfect_ics``, and rerunning ``perfect_model_obs``. This execution of ``perfect_model_obs`` will advance the model
state from the end of the first 24,000 steps to the end of an additional 24,000 steps and place the final state in
``perfect_restart``. The rest of the namelists in ``input.nml`` should remain unchanged.

::


   &perfect_model_obs_nml
      start_from_restart    = .true.,
      output_restart        = .true.,
      async                 = 0,
      init_time_days        = 0,
      init_time_seconds     = 0,
      first_obs_days        = -1,
      first_obs_seconds     = -1,
      last_obs_days         = -1,
      last_obs_seconds      = -1,
      output_interval       = 1,
      restart_in_file_name  = "perfect_ics",
      restart_out_file_name = "perfect_restart",
      obs_seq_in_file_name  = "obs_seq.in",
      obs_seq_out_file_name = "obs_seq.out",
      adv_ens_command       = "./advance_ens.csh"  /

::


   $ cp perfect_restart perfect_ics
   $ ./perfect_model_obs

A ``True_State.nc`` file is also created. It contains the 'true' state of the integration.

Generating the ensemble
^^^^^^^^^^^^^^^^^^^^^^^

This step (#2 from above) is done with the program ``filter``, which also uses the Fortran90 namelist mechanism for
input. It is now necessary to copy the ``input.nml.filter_default`` namelist to ``input.nml``.

::


   $ cp input.nml.filter_default
   $ input.nml

You may also build one master namelist containting all the required namelists. Having unused namelists in the
``input.nml`` does not hurt anything, and it has been so useful to be reminded of what is possible that we made it an
error to NOT have a required namelist. Take a peek at any of the other models for examples of a "fully qualified"
``input.nml``.

*HINT:* if you used ``svn`` to get the project, try 'svn revert input.nml' to restore the namelist that was distributed
with the project - which DOES have all the namelist blocks. Just be sure the values match the examples here.

.. container:: routineIndent1

   ::

      &filter_nml
         async                    = 0,
         adv_ens_command          = "./advance_model.csh",
         ens_size                 = 100,
         start_from_restart       = .false.,
         output_restart           = .true.,
         obs_sequence_in_name     = "obs_seq.out",
         obs_sequence_out_name    = "obs_seq.final",
         restart_in_file_name     = "perfect_ics",
         restart_out_file_name    = "filter_restart",
         init_time_days           = 0,
         init_time_seconds        = 0,
         first_obs_days           = -1,
         first_obs_seconds        = -1,
         last_obs_days            = -1,
         last_obs_seconds         = -1,
         num_output_state_members = 20,
         num_output_obs_members   = 20,
         output_interval          = 1,
         num_groups               = 1,
         input_qc_threshold       =  4.0,
         outlier_threshold        = -1.0,
         output_forward_op_errors = .false.,
         output_timestamps        = .false.,
         output_inflation         = .true.,

         inf_flavor               = 0,                       0,
         inf_start_from_restart   = .false.,                 .false.,
         inf_output_restart       = .false.,                 .false.,
         inf_deterministic        = .true.,                  .true.,
         inf_in_file_name         = 'not_initialized',       'not_initialized',
         inf_out_file_name        = 'not_initialized',       'not_initialized',
         inf_diag_file_name       = 'not_initialized',       'not_initialized',
         inf_initial              = 1.0,                     1.0,
         inf_sd_initial           = 0.0,                     0.0,
         inf_lower_bound          = 1.0,                     1.0,
         inf_upper_bound          = 1000000.0,               1000000.0,
         inf_sd_lower_bound       = 0.0,                     0.0
      /

      &smoother_nml
         num_lags              = 0,
         start_from_restart    = .false.,
         output_restart        = .false.,
         restart_in_file_name  = 'smoother_ics',
         restart_out_file_name = 'smoother_restart'  /

      &ensemble_manager_nml
         single_restart_file_in  = .true.,
         single_restart_file_out = .true.,
         perturbation_amplitude  = 0.2  /

      &assim_tools_nml
         filter_kind                     = 1,
         cutoff                          = 0.2,
         sort_obs_inc                    = .false.,
         spread_restoration              = .false.,
         sampling_error_correction       = .false.,
         adaptive_localization_threshold = -1,
         print_every_nth_obs             = 0  /

      &cov_cutoff_nml
         select_localization = 1  /

      &reg_factor_nml
         select_regression    = 1,
         input_reg_file       = "time_mean_reg",
         save_reg_diagnostics = .false.,
         reg_diagnostics_file = "reg_diagnostics"  /

      &obs_sequence_nml
         write_binary_obs_sequence = .false.  /

      &obs_kind_nml
         assimilate_these_obs_types = 'RAW_STATE_VARIABLE'  /

      &assim_model_nml
         write_binary_restart_files = .true. /

      &model_nml
         sigma  = 10.0,
         r      = 28.0,
         b      = 2.6666666666667,
         deltat = 0.01,
         time_step_days = 0,
         time_step_seconds = 3600  /

      &utilities_nml
         TERMLEVEL = 1,
         logfilename = 'dart_log.out'  /

Only the non-obvious(?) entries for ``filter_nml`` will be discussed.

+------------------------------+--------------------------------------------------------------------------------------+
| namelist variable            | description                                                                          |
+==============================+======================================================================================+
| ``ens_size``                 | Number of ensemble members. 100 is sufficient for most of the L63 exercises.         |
+------------------------------+--------------------------------------------------------------------------------------+
| ``start_from_restart``       | when '.false.', ``filter`` will generate its own ensemble of initial conditions. It  |
|                              | is important to note that the filter still makes use of the file named by            |
|                              | ``restart_in_file_name`` (i.e. ``perfect_ics``) by randomly perturbing these state   |
|                              | variables.                                                                           |
+------------------------------+--------------------------------------------------------------------------------------+
| ``num_output_state_members`` | specifies the number of state vectors contained in the netCDF diagnostic files. May  |
|                              | be a value from 0 to ``ens_size``.                                                   |
+------------------------------+--------------------------------------------------------------------------------------+
| ``num_output_obs_members``   | specifies the number of 'observations' (derived from applying the forward operator   |
|                              | to the state vector) are contained in the ``obs_seq.final`` file. May be a value     |
|                              | from 0 to ``ens_size``                                                               |
+------------------------------+--------------------------------------------------------------------------------------+
| ``inf_flavor``               | A value of 0 results in no inflation.(spin-up)                                       |
+------------------------------+--------------------------------------------------------------------------------------+

The filter is told to generate its own ensemble initial conditions since ``start_from_restart`` is '.false.'. However,
it is important to note that the filter still makes use of ``perfect_ics`` which is set to be the
``restart_in_file_name``. This is the model state generated from the first 24,000 step model integration by
``perfect_model_obs``. ``Filter`` generates its ensemble initial conditions by randomly perturbing the state variables
of this state.

``num_output_state_members`` are '.true.' so the state vector is output at every time for which there are observations
(once a day here). ``Posterior_Diag.nc`` and ``Prior_Diag.nc`` then contain values for 20 ensemble members once a day.
Once the namelist is set, execute ``filter`` to integrate the ensemble forward for 24,000 steps with the final ensemble
state written to the ``filter_restart``. Copy the ``perfect_model_obs`` restart file ``perfect_restart`` (the \`true
state') to ``perfect_ics``, and the ``filter`` restart file ``filter_restart`` to ``filter_ics`` so that future
assimilation experiments can be initialized from these spun-up states.

.. container:: unix

   ::

      ./filter
      cp perfect_restart perfect_ics
      cp filter_restart filter_ics

The spin-up of the ensemble can be viewed by examining the output in the netCDF files ``True_State.nc`` generated by
``perfect_model_obs`` and ``Posterior_Diag.nc`` and ``Prior_Diag.nc`` generated by ``filter``. To do this, see the
detailed discussion of matlab diagnostics in Appendix I.

3. Simulate a particular observing system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Begin by using ``create_obs_sequence`` to generate an observation set in which each of the 3 state variables of L63 is
observed with an observational error variance of 1.0 for each observation. To do this, use the following input sequence
(the text including and after # is a comment and does not need to be entered):

============= ===========================================================
*4*           # upper bound on num of observations in sequence
*0*           # number of copies of data (0 for just a definition)
*0*           # number of quality control values per field (0 or greater)
*0*           # -1 to exit/end observation definitions
*-1*          # observe state variable 1
*0 0*         # time -- days, seconds
*1.0*         # observational variance
*0*           # -1 to exit/end observation definitions
*-2*          # observe state variable 2
*0 0*         # time -- days, seconds
*1.0*         # observational variance
*0*           # -1 to exit/end observation definitions
*-3*          # observe state variable 3
*0 0*         # time -- days, seconds
*1.0*         # observational variance
*-1*          # -1 to exit/end observation definitions
*set_def.out* # Output file name
============= ===========================================================

Now, generate an observation sequence definition by running ``create_fixed_network_seq`` with the following input
sequence:

============= ===============================================================
*set_def.out* # Input observation set definition file
*1*           # Regular spaced observation interval in time
*1000*        # 1000 observation times
*0, 43200*    # First observation after 12 hours (0 days, 12 \* 3600 seconds)
*0, 43200*    # Observations every 12 hours
*obs_seq.in*  # Output file for observation sequence definition
============= ===============================================================

4. Generate a particular observing system and true state
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An observation sequence file is now generated by running ``perfect_model_obs`` with the namelist values (unchanged from
step 2):

.. container:: routineIndent1

   ::

      &perfect_model_obs_nml
         start_from_restart    = .true.,
         output_restart        = .true.,
         async                 = 0,
         init_time_days        = 0,
         init_time_seconds     = 0,
         first_obs_days        = -1,
         first_obs_seconds     = -1,
         last_obs_days         = -1,
         last_obs_seconds      = -1,
         output_interval       = 1,
         restart_in_file_name  = "perfect_ics",
         restart_out_file_name = "perfect_restart",
         obs_seq_in_file_name  = "obs_seq.in",
         obs_seq_out_file_name = "obs_seq.out",
         adv_ens_command       = "./advance_ens.csh"  /

This integrates the model starting from the state in ``perfect_ics`` for 1000 12-hour intervals outputting synthetic
observations of the three state variables every 12 hours and producing a netCDF diagnostic file, ``True_State.nc``.

5. Filtering
~~~~~~~~~~~~

Finally, ``filter`` can be run with its namelist set to:

.. container:: routineIndent1

   ::

      &filter_nml
         async                    = 0,
         adv_ens_command          = "./advance_model.csh",
         ens_size                 = 100,
         start_from_restart       = .true.,
         output_restart           = .true.,
         obs_sequence_in_name     = "obs_seq.out",
         obs_sequence_out_name    = "obs_seq.final",
         restart_in_file_name     = "filter_ics",
         restart_out_file_name    = "filter_restart",
         init_time_days           = 0,
         init_time_seconds        = 0,
         first_obs_days           = -1,
         first_obs_seconds        = -1,
         last_obs_days            = -1,
         last_obs_seconds         = -1,
         num_output_state_members = 20,
         num_output_obs_members   = 20,
         output_interval          = 1,
         num_groups               = 1,
         input_qc_threshold       =  4.0,
         outlier_threshold        = -1.0,
         output_forward_op_errors = .false.,
         output_timestamps        = .false.,
         output_inflation         = .true.,

         inf_flavor               = 0,                       0,
         inf_start_from_restart   = .false.,                 .false.,
         inf_output_restart       = .false.,                 .false.,
         inf_deterministic        = .true.,                  .true.,
         inf_in_file_name         = 'not_initialized',       'not_initialized',
         inf_out_file_name        = 'not_initialized',       'not_initialized',
         inf_diag_file_name       = 'not_initialized',       'not_initialized',
         inf_initial              = 1.0,                     1.0,
         inf_sd_initial           = 0.0,                     0.0,
         inf_lower_bound          = 1.0,                     1.0,
         inf_upper_bound          = 1000000.0,               1000000.0,
         inf_sd_lower_bound       = 0.0,                     0.0
       /

``filter`` produces two output diagnostic files, ``Prior_Diag.nc`` which contains values of the ensemble mean, ensemble
spread, and ensemble members for 12- hour lead forecasts before assimilation is applied and ``Posterior_Diag.nc`` which
contains similar data for after the assimilation is applied (sometimes referred to as analysis values).

Now try applying all of the matlab diagnostic functions described in the Matlab® Diagnostics section.

The tutorial
------------

The ``DART/tutorial`` documents are an excellent way to kick the tires on DART and learn about ensemble data
assimilation. If you have gotten this far, you can run anything in the tutorial.

Matlab® diagnostics
-------------------

The output files are netCDF files, and may be examined with many different software packages. We happen to use Matlab®,
and provide our diagnostic scripts in the hopes that they are useful.

The diagnostic scripts and underlying functions reside in two places: ``DART/diagnostics/matlab`` and ``DART/matlab``.
They are reliant on the public-domain `netcdf
toolbox <http://woodshole.er.usgs.gov/staffpages/cdenham/public_html/MexCDF/nc4ml5.html>`__ from
``http://woodshole.er.usgs.gov/staffpages/cdenham/public_html/MexCDF/nc4ml5.html`` as well as the public-domain `CSIRO
matlab/netCDF interface <http://www.marine.csiro.au/sw/matlab-netcdf.html>`__ from
``http://www.marine.csiro.au/sw/matlab-netcdf.html``. If you do not have them installed on your system and want to use
Matlab to peruse netCDF, you must follow their installation instructions. The 'interested reader' may want to look at
the ``DART/matlab/startup.m`` file I use on my system. If you put it in your ``$HOME/matlab`` directory, it is invoked
every time you start up Matlab.

| Once you can access the ``getnc`` function from within Matlab, you can use our diagnostic scripts. It is necessary to
  prepend the location of the ``DART/matlab`` scripts to the ``matlabpath``. Keep in mind the location of the netcdf
  operators on your system WILL be different from ours ... and that's OK.

.. container:: unix

   ::

      [models/lorenz_63/work]$ matlab -nojvm

                                                   < M A T L A B >
                                       Copyright 1984-2002 The MathWorks, Inc.
                                           Version 6.5.0.180913a Release 13
                                                     Jun 18 2002

        Using Toolbox Path Cache.  Type "help toolbox_path_cache" for more info.
       
        To get started, type one of these: helpwin, helpdesk, or demo.
        For product information, visit www.mathworks.com.

      >> which getnc
      /contrib/matlab/matlab_netcdf_5_0/getnc.m
      >>ls *.nc

      ans =

      Posterior_Diag.nc  Prior_Diag.nc  True_State.nc


      >>path('../../matlab',path)
      >>path('../../diagnostics/matlab',path)
      >>which plot_ens_err_spread
      ../../matlab/plot_ens_err_spread.m
      >>help plot_ens_err_spread

        DART : Plots summary plots of the ensemble error and ensemble spread.
                               Interactively queries for the needed information.
                               Since different models potentially need different 
                               pieces of information ... the model types are 
                               determined and additional user input may be queried.
       
        Ultimately, plot_ens_err_spread will be replaced by a GUI.
        All the heavy lifting is done by PlotEnsErrSpread.
       
        Example 1 (for low-order models)
       
        truth_file = 'True_State.nc';
        diagn_file = 'Prior_Diag.nc';
        plot_ens_err_spread

      >>plot_ens_err_spread

And the matlab graphics window will display the spread of the ensemble error for each state variable. The scripts are
designed to do the "obvious" thing for the low-order models and will prompt for additional information if needed. The
philosophy of these is that anything that starts with a lower-case *plot\_\ some_specific_task* is intended to be
user-callable and should handle any of the models. All the other routines in ``DART/matlab`` are called BY the
high-level routines.

+-------------------------------+-------------------------------------------------------------------------------------+
| Matlab script                 | description                                                                         |
+===============================+=====================================================================================+
| ``plot_bins``                 | plots ensemble rank histograms                                                      |
+-------------------------------+-------------------------------------------------------------------------------------+
| ``plot_correl``               | Plots space-time series of correlation between a given variable at a given time and |
|                               | other variables at all times in a n ensemble time sequence.                         |
+-------------------------------+-------------------------------------------------------------------------------------+
| ``plot_ens_err_spread``       | Plots summary plots of the ensemble error and ensemble spread. Interactively        |
|                               | queries for the needed information. Since different models potentially need         |
|                               | different pieces of information ... the model types are determined and additional   |
|                               | user input may be queried.                                                          |
+-------------------------------+-------------------------------------------------------------------------------------+
| ``plot_ens_mean_time_series`` | Queries for the state variables to plot.                                            |
+-------------------------------+-------------------------------------------------------------------------------------+
| ``plot_ens_time_series``      | Queries for the state variables to plot.                                            |
+-------------------------------+-------------------------------------------------------------------------------------+
| ``plot_phase_space``          | Plots a 3D trajectory of (3 state variables of) a single ensemble member.           |
|                               | Additional trajectories may be superimposed.                                        |
+-------------------------------+-------------------------------------------------------------------------------------+
| ``plot_total_err``            | Summary plots of global error and spread.                                           |
+-------------------------------+-------------------------------------------------------------------------------------+
| ``plot_var_var_correl``       | Plots time series of correlation between a given variable at a given time and       |
|                               | another variable at all times in an ensemble time sequence.                         |
+-------------------------------+-------------------------------------------------------------------------------------+

Bias, filter divergence and covariance inflation (with the l63 model)
---------------------------------------------------------------------

One of the common problems with ensemble filters is filter divergence, which can also be an issue with a variety of
other flavors of filters including the classical Kalman filter. In filter divergence, the prior estimate of the model
state becomes too confident, either by chance or because of errors in the forecast model, the observational error
characteristics, or approximations in the filter itself. If the filter is inappropriately confident that its prior
estimate is correct, it will then tend to give less weight to observations than they should be given. The result can be
enhanced overconfidence in the model's state estimate. In severe cases, this can spiral out of control and the ensemble
can wander entirely away from the truth, confident that it is correct in its estimate. In less severe cases, the
ensemble estimates may not diverge entirely from the truth but may still be too confident in their estimate. The result
is that the truth ends up being farther away from the filter estimates than the spread of the filter ensemble would
estimate. This type of behavior is commonly detected using rank histograms (also known as Talagrand diagrams). You can
see the rank histograms for the L63 initial assimilation by using the matlab script ``plot_bins``.

A simple, but surprisingly effective way of dealing with filter divergence is known as covariance inflation. In this
method, the prior ensemble estimate of the state is expanded around its mean by a constant factor, effectively
increasing the prior estimate of uncertainty while leaving the prior mean estimate unchanged. The program ``filter`` has
a group of namelist parameters that controls the application of covariance inflation. For a simple set of inflation
values, you will set ``inf_flavor``, and ``inf_initial``. These values come in pairs; the first value controls inflation
of the prior ensemble values, while the second controls inflation of the posterior values. Up to this point
``inf_flavor`` has been set to 0 indicating that the prior ensemble is left unchanged. Setting the first value of
``inf_flavor`` to 3 enables one variety of inflation. Set ``inf_initial`` to different values (try 1.05 and 1.10 and
other values). In each case, use the diagnostic matlab tools to examine the resulting changes to the error, the ensemble
spread (via rank histogram bins, too), etc. What kind of relation between spread and error is seen in this model?

There are many more options for inflation, including spatially and temporarily varying values, with and without damping.
See the discussion of all inflation-related namelist items
`Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/filter/filter.html#Inflation>`__ or `local
file <../../filter/filter.html#Inflation>`__.

Synthetic observations
----------------------

Synthetic observations are generated from a \`perfect' model integration, which is often referred to as the \`truth' or
a \`nature run'. A model is integrated forward from some set of initial conditions and observations are generated as *y
= H(x) + e* where *H* is an operator on the model state vector, *x*, that gives the expected value of a set of
observations, *y*, and *e* is a random variable with a distribution describing the error characteristics of the
observing instrument(s) being simulated. Using synthetic observations in this way allows students to learn about
assimilation algorithms while being isolated from the additional (extreme) complexity associated with model error and
unknown observational error characteristics. In other words, for the real-world assimilation problem, the model has
(often substantial) differences from what happens in the real system and the observational error distribution may be
very complicated and is certainly not well known. Be careful to keep these issues in mind while exploring the
capabilities of the ensemble filters with synthetic observations.

Notes for current users
-----------------------

If you have been updating from the head of the DART subversion repository (the "trunk") you will not notice much
difference between that and the Kodiak release. If you are still running the Jamaica release there are many new models,
new observation types, capabilities in the assimilation tools, new diagnostics, and new utilities. There is a very short
list of non-backwards compatible changes (see below), and then a long list of new options and functions.

In recent years we have been adding new functionality to the head of the subversion trunk and just testing it and
keeping it in working order, maintaining backwards compatibility. We now have many development tasks which will require
non-compatible changes in interfaces and behavior. Further DART development will occur on a branch, so checking out
either the Kodiak branch or the head of the repository is the recommended way to update your DART tree.

Non-backwards compatible changes
--------------------------------

Changes in the Kodiak release which are *not* backwards compatible with the Jamaica release (svn revision number 2857,
12 April 2007):

#. &filter_nml used to have a single entry to control whether to read in both the inflation values and standard
   deviations from a file or use the settings in the namelist. The old namelist item, ``inf_start_from_file``, has been
   replaced by two items that allow the inflation values and the standard deviation to be read in separately. The new
   namelist items are ``inf_initial_from_file`` and ``inf_sd_initial_from_file``. See the filter namelist documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/filter/filter.html#Namelist>`__ or `local
   file <../../filter/filter.html#Namelist>`__ for more details.

#. The WRF/DART converter program used to be called ``dart_tf_wrf``, had no namelist, and you entered ``T`` or ``F`` to
   indicate which direction you were converting. Now we have ``dart_to_wrf`` and ``wrf_to_dart`` (documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/wrf/models/wrf/WRF_DART_utilities/dart_to_wrf.html>`__)
   each with a namelist to control various options.

#. The CAM/DART converter programs used to be called ``trans_sv_pv`` and ``trans_pv_sv``, with no namelists, and with
   several specialized variants (e.g. ``trans_pv_sv_time0``). Now we have ``cam_to_dart`` (documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/cam/cam_to_dart.html>`__ ) and ``dart_to_cam`` (documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/cam/dart_to_cam.html>`__ ) each with a namelist to control various options.

#. The ``obs_def_radar_mod.f90`` radar observation module was completely rewritten and the namelist has changed
   substantially. See the module documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/forward_operators/obs_def_radar_mod.html>`__
   or :doc:`../../observations/forward_operators/obs_def_radar_mod` for details. For example, the defaults for the
   old code were:

   ::

      &obs_def_radar_mod_nml
         convert_to_dbz            =  .true. ,
         dbz_threshold             =   0.001 ,
         apply_ref_limit_to_obs    = .false. ,
         reflectivity_limit_obs    =     0.0 ,
         lowest_reflectivity_obs   = -888888.0,
         apply_ref_limit_to_state  = .false. ,
         reflectivity_limit_state  =     0.0 ,
         lowest_reflectivity_state = -888888.0 /

   and the new ones are:

   ::

      &obs_def_radar_mod_nml
         apply_ref_limit_to_obs     =  .true. ,
         reflectivity_limit_obs     =     0.0 ,
         lowest_reflectivity_obs    =     0.0 ,
         apply_ref_limit_to_fwd_op  =  .true. ,
         reflectivity_limit_fwd_op  =     0.0 ,
         lowest_reflectivity_fwd_op =     0.0 ,
         dielectric_factor          =   0.224 ,
         n0_rain                    =   8.0e6 ,
         n0_graupel                 =   4.0e6 ,
         n0_snow                    =   3.0e6 ,
         rho_rain                   =  1000.0 ,
         rho_graupel                =   400.0 ,
         rho_snow                   =   100.0 ,
         allow_wet_graupel          = .false.,
         microphysics_type          =       3 ,
         allow_dbztowt_conv         = .false. /

#. The WRF &model_mod namelist has changed. It now requires a ``wrf_state_variables`` list to choose which WRF fields
   are put into the state vector. The order of the field names in the list sets the order of the fields in the state
   vector. See the WRF model_mod documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/wrf/model_mod.html#Namelist>`__ or `local
   file <../../models/wrf/model_mod.html#Namelist>`__ for details. Although they haven't been removed from the
   namelist, the following items have no effect on the code anymore:

   -  num_moist_vars
   -  surf_obs
   -  soil_data
   -  h_diab

#. The WRF model_mod now computes geometric heights instead of geopotential heights. It also uses the staggered grids as
   read in from the ``wrfinput_dNN`` file(s) instead of interpolating in the non-staggered grid to get individual cell
   corners.

#. The code in ``filter.f90`` was corrected to match the documentation for how the namelist item ``input_qc_threshold``
   is handled. (See filter namelist documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/filter/filter.html#Namelist>`__ or `local
   file <../../filter/filter.html#Namelist>`__.) In the Jamaica release, observations with incoming data QC values
   greater than or equal to the namelist setting were discarded. Now only incoming data QC values greater than the
   ``input_qc_threshold`` are discarded (values equal to the threshold are now kept).

#. The ``merge_obs_seq`` utility has been replaced by the more comprehensive ``obs_sequence_tool`` utility. See the
   documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/programs/obs_sequence_tool/assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html>`__
   or
   :doc:`../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`.

#. The prepbufr observation converter was located in the ``DART/ncep_obs`` directory in the last release. It has been
   moved to be with the other programs that convert various types of observation files into DART format. It is now
   located in ``DART/observations/NCEP``.

#. The sampling error correction generator program in ``DART/system_simulation`` now has a namelist &full_error_nml. See
   the documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/programs/system_simulation/system_simulation.html>`__
   or :doc:`../../assimilation_code/programs/system_simulation/system_simulation` for more details. Tables for 40
   common ensemble sizes are pregenerated in the ``DART/system_simulation/final_full_precomputed_tables`` directory, and
   instructions for generating tables for other ensemble sizes are given.

#. Most ``work`` directories now have a ``quickbuild.csh`` script which recompiles all the executables instead of a
   ``workshop_setup.csh`` script. (Those directories used in the tutorial have both.) To control whether ``filter`` is
   compiled with or without MPI (as a parallel program or not) the ``quickbuild.csh`` script takes the optional
   arguments ``-mpi`` or ``-nompi``.

#. The ``preprocess`` program was changed so that any obs_def files with module definitions are directly included in the
   single ``obs_def_mod.f90`` file. This means that as you add and delete obs_def modules from your &preprocess_nml
   namelist and rerun ``preprocess`` you no longer have to add and delete different obs_def modules from your
   ``path_names_*`` files.

#. The utilities module now calls a function in the mpi_utilities code to exit MPI jobs cleanly. This requires that
   non-mpi programs now include the ``null_mpi_utilities_mod.f90`` file in their ``path_names_*`` files.

#. The ``DART/mpi_utilities`` directory as distributed now works with all compilers except for gfortran. In
   ``DART/mpi_utilities`` is a ``./fixsystem`` script that when executed will change the source files so they will
   compile with gfortran. Previous releases compiled with gfortran as distributed but no other compilers.

#. The GPS Radio Occultation observation forward operator code now requires a namelist, ``&obs_def_gps_nml``. See the
   GPS documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/forward_operators/obs_def_gps_mod.html#Namelist>`__
   or `local file <../../observations/forward_operators/obs_def_gps_mod.html#Namelist>`__ for details on what to add.
   All ``input.nml`` files in the repository have had this added if they have the GPS module in their
   ``&preprocess_nml`` namelist.

New features
------------

-  Inflation Damping

   -  Handles the case where observation density is irregular in time, e.g. areas which were densely observed at one
      point are no longer observed. Adaptive inflation values can grow large where the observations are dense, and if
      that region is no longer observed the inflation is not recomputed. Inflation damping shrinks the inflation values
      and compensates for this. See the inflation documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/filter/filter.html#Inflation>`__ or `local
      file <../../filter/filter.html#Inflation>`__ for more details and paper references.

-  Sampling Error Correction

   -  Compensates for the numbers of ensembles being small compared to the number of degrees of freedom in the system.
      See the last item in this section of the documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/filter/filter.html#GettingStarted>`__ or `local
      file <../../filter/filter.html#GettingStarted>`__ for more details.

-  Adaptive Localization and Localization Diagnostics

   -  See a discussion of localization-related issues
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assim_tools/assim_tools_mod.html#Localization>`__
      or `local file <../../assim_tools/assim_tools_mod.html#Localization>`__.

-  Scale height vertical localization option in 3d models

   -  See a discussion of specifying vertical localization in terms of scale height
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/location/threed_sphere/location_mod.html#Namelist>`__
      or `local file <../../location/threed_sphere/location_mod.html#Namelist>`__. See the `Wikipedia
      page <http://en.wikipedia.org/wiki/Scale_height>`__ for a discussion of how scale height is defined. Note that
      there is no support in the diagnostic Matlab routines for observations using scale height as the vertical
      coordinate.

-  CAM supports FV code, PBS scripting

   -  See details on the features of the CAM/DART system
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/cam/model_mod.html>`__ .

-  Boxcar Kernel Filter Option

   -  See how to select this filter option in the namelist
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assim_tools/assim_tools_mod.html#FilterTypes>`__
      or `local file <../../assim_tools/assim_tools_mod.html#FilterTypes>`__.

-  Option for "undefined vertical location" for obs using the 3d sphere locations

   -  See how to specify this option when creating observations
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/location/threed_sphere/location_mod.html>`__ or
      :doc:`../../assimilation_code/location/threed_sphere/location_mod`.

-  Schedule module for repeated time intervals

   -  See documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/modules/utilities/schedule_mod.html>`__
      or :doc:`../../assimilation_code/modules/utilities/schedule_mod`.

-  Support for 2 different Mars calendars in time manager

   -  Gregorian Mars
   -  Solar Mars

-  Code corrections to make the smoother run correctly
-  Forward operators now have access to the ensemble number and the state time if they want to make use of this
   information

   -  See the "Get Expected Obs From Def" section of the obs_def documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/forward_operators/observations/forward_operators/obs_def_mod.html>`__
      for details on how to use these values. This change is fully backwards-compatible with existing forward operator code.

-  Option to output all echo of namelist values to a separate log file

   -  See the utilities module documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/modules/utilities/utilities_mod.html#Namelist>`__
      or `local file <../../assimilation_code/modules/utilities/utilities_mod.html#Namelist>`__ for how to select
      where the contents of all namelists are output.

-  Large file support for netCDF

   -  See the `Unidata netCDF documentation <http://www.unidata.ucar.edu/software/netcdf/faq-lfs.html>`__ pages for more
      information about what large file support gives you and what it is compatible with.

-  Better support for adaptive localization

   -  See the Localization section of the assim_tools documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assim_tools/assim_tools_mod.html#Localization>`__
      or `local file <../../assim_tools/assim_tools_mod.html#Localization>`__ for details.

-  Option to localize with different distances based on observation type

   -  See the Localization section of the assim_tools documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assim_tools/assim_tools_mod.html#Localization>`__
      or `local file <../../assim_tools/assim_tools_mod.html#Localization>`__ for details.

-  The error handler can take up to 3 lines of text so you can give more informative error messages on exit

   -  See the utilities module documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/modules/utilities/utilities_mod.html#Interface>`__
      or `local file <../../assimilation_code/modules/utilities/utilities_mod.html#Interface>`__ for details.

-  Option to output ensemble mean in restart file format when filter exits

   -  See the filter program namelist documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/filter/filter.html#Namelist>`__ or `local
      file <../../filter/filter.html#Namelist>`__ for details.

-  The start of a suite of forecast verification and evaluation tools

   -  See the verification tool documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/obs_sequence/assimilation_code/programs/obs_seq_verify/obs_seq_verify.html>`__
      or :doc:`../../assimilation_code/programs/obs_seq_verify/obs_seq_verify` for details.

-  Performance improvement in the internal transposes for very large state vectors. all_vars_to_all_copies() now has a
   single receiver and multiple senders, which is much faster than the converse.
-  Better support for users who redefine R8 to be R4, so that filter runs in single precision. Fixed code which was
   technically correct but numerically unstable in single precision when computing variance and covariances.
-  Fixed a case in the 3D sphere locations code which made it possible that some observations and state variables at
   higher latitudes might not be impacted by observations which were barely within the localization cutoff.
-  The observation type table at the top of all obs_seq files now only contains the types actually found in the file.
-  When one or more ensemble members fail to compute a valid forward operator, the prior and/or posterior mean and
   standard deviation will be set to MISSING_R8 in the output obs_seq.final file in addition to setting the DART QC
   flag.
-  Use less stack space by allocating large arrays instead of declaring them as local (stack) variables in routines
-  The copyright has changed from GPL (GNU) to an NCAR-specific one which is found
   `here <http://www.image.ucar.edu/DAReS/DART/DART_download>`__.

New models
----------

-  POP Ocean Model

   -  DART interface documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/POP/model_mod.html>`__ or
      :doc:`../../models/POP/readme`. Documentation for the model itself `in
      CESM <http://www.cesm.ucar.edu/models/ccsm2.0.1/pop/>`__ and `stand-alone version from Los
      Alamos <http://climate.lanl.gov/Models/POP/>`__.

-  NCOMMAS Mesoscale Atmospheric Model

   -  DART interface documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/NCOMMAS/model_mod.html>`__ or
      :doc:`../../models/NCOMMAS/readme`. Documentation for the model itself from NSSL, Norman, OK. is at
      `NCOMMAS <http://code.google.com/p/enkf-nssl-commas>`__.

-  COAMPS Atmosphere Model

   -  Dart interface documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/coamps_next/model_mod.html>`__ or
      :doc:`../../models/coamps_nest/readme`. Documentation for the model itself is at
      `COAMPS <http://www.nrlmry.navy.mil/coamps-web/web/home>`__. The original version of the COAMPS interface code and
      scripts was contributed by Tim Whitcomb, NRL, Monterey. An updated version was contributed by Alex Reinecke, NRL,
      Monterey.
      The primary differences between the current version and the original COAMPS model code are:

      -  the ability to assimilate nested domains
      -  assimilates real observations
      -  a simplified way to specify the state vector
      -  I/O COAMPS data files
      -  extensive script updates to accommodate additional HPC environments

-  NOGAPS Global Atmosphere Model

   -  The Navy's operational global atmospheric prediction system. See
      `here <http://www.srh.noaa.gov/ssd/nwpmodel/html/nogover.htm>`__ for an overview of the 4.0 version of the model.
      For more information on the NOGAPS/DART system, contact Jim Hansen, jim.hansen at nrlmry.navy.mil

-  AM2 Atmosphere Model

   -  Dart interface documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/AM2/model_mod.html>`__ or
      :doc:`../../models/am2/readme`. The GFDL atmosphere model documentation is at
      `AM2 <http://data1.gfdl.noaa.gov/~arl/pubrel/m/am2/doc/>`__.

-  MIT Global Ocean Model

   -  Dart interface documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/MITgcm_ocean/model_mod.html>`__ or
      :doc:`../../models/MITgcm_ocean/readme`. The `ocean
      component <http://paoc2001.mit.edu/cmi/development/ocean.htm>`__ of the MIT global model suite.

-  Simple Advection Model

   -  Dart interface documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/simple_advection/model_mod.html>`__ or
      :doc:`../../models/simple_advection/readme`. A simple model of advecting tracers such as CO.

-  Global/Planet WRF

   -  A version of the WRF weather model adapted for `global use or for other planets <http://planetwrf.com/>`__.

-  TIEgcm Thermosphere/Ionosphere Model

   -  Dart interface documentation
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/TIEgcm/model_mod.html>`__ or
      :doc:`../../models/tiegcm/readme`. Documentation for the thermosphere and ionosphere model from the NCAR HAO
      (High Altitude Observatory) Division is at `TIEgcm <http://cism.hao.ucar.edu/models_tiegcm.html>`__.

The ``DART/models/template`` directory contains sample files for adding a new model. See `this
section <http://www.image.ucar.edu/DAReS/DART/DART_Documentation.php#adding_a_model>`__ of the DART web pages for more
help on adding a new model.

Changed models
--------------

-  WRF

   -  The WRF fields in the DART state vector are namelist settable, with the order of the names in the namelist
      controlling the order in the state vector. No assumptions are made about number of moist variables; all WRF fields
      must be listed explicitly. The conversion tools dart_to_wrf and wrf_to_dart (Documented here
      `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/wrf/models/wrf/WRF_DART_utilities/dart_to_wrf.html>`__ )
      use this same namelist, so it is simpler to avoid mismatches between the DART restart files and what the WRF model_mod is expecting.
   -  Support for the single column version of WRF has been incorporated into the standard WRF model_mod.
   -  advance_model.csh script reworked by Josh Hacker, Ryan Torn, and Glen Romine to add function and simplify the
      script. It now supports a restart-file-per-member, simplifies the time computations by using the advance_time
      executable, and handles boundary files more cleanly. Plus added many additional comments, and ways to select
      various options by setting shell variables at the top of the script.
   -  Updates from Tim and Glen:
      - Changed the variable name for the longitude array to better match that used in WRF: XLON_d0\* to XLONG_d0\*
      - Added the staggered coordinate variables (XLONG_U_d0*, XLAT_U_d0*, XLONG_V_d0*, XLAT_V_d0*, ZNW_d0*)
      - Use the staggered variables to look up point locations when interpolating in a staggered grid. Old code used to
      compute the staggered points from the unstaggered grid, which was slightly inaccurate.
      - Added additional attribute information, supporting long_name, description (same info as long_name which is the
      standard, but WRF calls this attribute 'description'), units (previously supported) and named coordinates for the
      X and Y directions (in keeping with WRF, we do not name the vertical coordinate).
   -  New scripts to generate LBC (lateral boundary condition) files for WRF runs.
-  CAM

   -  support for versions 4 and 5 of CAM, including tar files of changes that must be made to the CAM source tree and
      incorporated into the CAM executable
   -  support leap years
   -  use CLM restart file instead of initial file
   -  various scripting changes to support archiving
   -  save information from CAM for ocean and land forcing
   -  scripts to archive months of obs_seq.finals
   -  Added the changes needed to the CAM build tree for CAM 4.0.x
   -  Updates to CAM documentation to bring it in sync with the current code.
   -  All trans routines replaced with: dart_to_cam, cam_to_dart, and advance_time.
   -  Minor changes to CAM model_mod, including adding a routine to write out the times file so utilities can call it in
      a single location, plus additional optional arg on write routine.
   -  Most debugging output is off by default; a new namelist item 'print_details' will re-enable the original output.
   -  Added build support for new tools (closest member, common subset, fill inflation) and removed for obsolete (merge
      obs).
   -  The original 'trans' build files and src are now in a 'deprecated' directory so if users need them for backwards
      compatibility, they are still available.
   -  The archive scripts are updated for the HPSS (hsi) and the MSS versions (msrcp) are removed.
   -  The shell_scripts and full_experiment scripts are updated.
-  Lorenz 2004/2005

   -  Fixed a bug in the model advance code which was doing an extra divide by 2, causing incorrect results.

New observation types/sources
-----------------------------

-  MADIS
   Converters for METAR, Mesonet, Rawinsondes, ACARS, Marine, and Satellite Wind observations. Optionally output
   moisture obs as specific humidity, relative humidity, and/or dewpoint obs. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/observations/obs_converters/MADIS/MADIS.html>`__ .
-  SSEC
   Convert Satellite Wind obs to DART format. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/observations/obs_converters/SSEC/SSEC.html>`__ .
-  AIRS
   Satellite observed Temperature and Moisture obs. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/observations/obs_converters/AIRS/AIRS.html>`__ .
-  QUIKscat
   Satellite observed surface winds. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/quikscat/quikscat.html>`__ .
-  GTSPP
   Ocean obs. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/observations/obs_converters/GTSPP/GTSPP.html>`__ .
-  WOD
   World Ocean Database (currently 2009) Temperature and Salinity obs. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/observations/obs_converters/WOD/WOD.html>`__ .
-  CODAR
   High frequency radar obs of ocean surface velocity. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/obs_def/obs_def_ocean_mod.f90>`__ or `local
   file <../../obs_def/obs_def_ocean_mod.f90>`__.
-  VAR
   Little-r and radar obs. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/obs_converters/var/var.html>`__ .
-  Text
   Reads text data files, a good template for converting obs stored in files without some kind of data library format
   (netCDF, HDF, etc). Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/obs_converters/observations/obs_converters/text/text_to_obs.html>`__ .
-  Altimeter
   Altimeter observation type available from a variety of sources. Forward operator code
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/obs_def/obs_def_altimeter_mod.f90>`__ or `local
   file <../../obs_def/obs_def_altimeter_mod.f90>`__.
-  Dewpoint
   Dewpoint observation type available from a variety of sources. Forward operator code
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/obs_def/obs_def_dewpoint_mod.f90>`__ or `local
   file <../../obs_def/obs_def_dewpoint_mod.f90>`__.
-  Dropsonde
   Dropsonde observation type available to allow these observations to be distinguished from standard Radiosondes. Type
   defined in code
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/obs_def/obs_def_reanalysis_bufr_mod.f90>`__ or
   `local file <../../obs_def/obs_def_reanalysis_bufr_mod.f90>`__.
-  TES Radiances
   TES satellite radiance observations of Mars. Forward operator code
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/obs_def/obs_def_TES_nadir_mod.f90>`__ or `local
   file <../../obs_def/obs_def_TES_nadir_mod.f90>`__.
-  Hurricane/Tropical Storm Vortex Position
   Storm location, minimum central pressure, and maximum windspeed. Currently only implemented in the WRF model_mod
   interface code. Code `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/wrf/model_mod.html>`__
   or :doc:`../../models/wrf/readme`.

All the observation converters have moved to their own top level directory ``observations``. See the overview
documentation
`Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/obs_converters/observations.html>`__ for
general information on creating observation files for use in the ensemble assimilation system.

The GPS forward operators aren't new with this release, but the code has been revised several times. In particular,
there is a new namelist to set the maximum number of GPS obs supported in a single execution of filter or the obs diag
program. Generally the default value is large enough for anything less than a couple days, but if you are running a
month or longer of diagnostics for a time series you can easily exceed the compiled in maximum. See the GPS
documentation for creating GPS observation files
`Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/GPS/gps.html>`__ or
:doc:`../../observations/obs_converters/gps/gps`, and the forward operator documentation
`Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/forward_operators/obs_def_gps_mod.html>`__
or :doc:`../../observations/forward_operators/obs_def_gps_mod`. There are also heavily revised scripts which download
and convert multiple days of GPS obs at a time, with options to delete downloaded files automatically. The scripts are
able to download GPS RO observations from any of about seven available satellites (in addition to the COSMIC array) from
the CDAAC web site.

There are two modules to set observation error values when creating new observation sequence files. One contains the
default values used by NCEP, and the other contains the values used by ECMWF. See the README file
`Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/observations/obs_error/README>`__ or `local
file <../../observations/obs_error/README>`__ for more details.

The radar module was completely overhauled and the namelist changed substantially. See the item above in the
non-backward compatible changes section for details.

The scripting for converting NCEP prepbufr observations has been improved. There are options to enable or disable the
'blocking' conversion, to create 6 hour or daily output files, to swap bytes for little-endian machines, and to run up
to a month of conversions in parallel if you have parallel hardware available.

There is a ``DART/observations/utilities`` directory where generic utilities can be built which are not dependent on any
particular model.

New diagnostics and documentation
---------------------------------

**Better Web Pages.** We've put a lot of effort into expanding our documentation. For example, please check out `the
Matlab diagnostics section <http://www.image.ucar.edu/DAReS/DART/DART_Documentation.php#mat_obs>`__ or the pages
outlining the `observation sequence file
contents <http://www.image.ucar.edu/DAReS/DART/DART_Observations.php#obs_seq_overview>`__.

But there's always more to add. **Please let us know where we are lacking.**

Other new stuff:

-  There is now a main ``index.html`` file
   (`Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/index.html>`__) in
   the DART distribution to quickly guide you to any of the documentation for the routines or modules.
-  DART_LAB
   Interactive Matlab GUI experiments and Powerpoint presentation of fundamental assimilation concepts
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/docs/DART_LAB/DART_LAB.html>`__ or
   :doc:`../DART_LAB/DART_LAB`.
-  link_obs.m
   Allows one to view multiple observation attributes simultaneously and dynamically select subsets of observations in
   one view and have those same obs highlighted in the other views. Commonly called 'data brushing'. Matlab source
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/link_obs.m>`__ or `local
   file <../../diagnostics/matlab/link_obs.m>`__.
-  obs_diag
   The ``obs_diag`` program has undergone extensive revision. User-defined levels for all coordinate
   (height/pressure/etc), arbitrary number of regions, the inclusion of separate copies for all DART QC values, can
   creates rank histograms from the ``obs_seq.final`` files, if possible, and more. See the documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/programs/obs_diag/oned/obs_diag.html%20assimilation_code/programs/obs_diag/threed_cartesian/obs_diag.html%20assimilation_code/programs/obs_diag/threed_sphere/obs_diag.html>`__ .
-  Comparing two (or more) experiments
   Matlab scripts to compare **multiple** (not just two) ``obs_diag_output.nc`` files on the same graphic to allow for
   easy examination of experiment attributes (rmse, biases, etc.). Some new utilities for subsetting observation
   sequence files in order to make fair comparisons are described below. Matlab source for ``two_experiments_profile.m``
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/two_experiments_profile.m>`__
   or `local file <../../diagnostics/matlab/two_experiments_profile.m>`__ and ``two_experiments_evolution.m``
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/two_experiments_evolution.m>`__
   or `local file <../../diagnostics/matlab/two_experiments_evolution.m>`__.
-  netCDF and Matlab
   The DART Matlab routines no longer depend on 4 third-party toolboxes; we are down to just
   `mexnc <http://mexcdf.sourceforge.net/downloads/>`__ and `snctools <http://mexcdf.sourceforge.net/downloads/>`__.
   Soon, we may just use snctools! See the documentation for how to configure Matlab to run the DART-provided scripts
   `Website <http://www.image.ucar.edu/DAReS/DART/DART_Documentation.php#configure_matlab>`__ or `local
   file <http://www.image.ucar.edu/DAReS/DART/DART_Documentation.php#configure_matlab>`__.
-  Matlab support for CAM.
   CAM is now fully supported for all the Matlab interfaces that are used in the demos - this includes the state-space
   tools in ``DART/matlab`` that allow for determining correlations among state variables, among other things.
-  Matlab support for WRF.
   WRF is now fully supported for all the Matlab interfaces that are used in the demos. This predominantly includes the
   state-space tools in the ``DART/matlab`` directory like ``plot_total_err``. The ``map_wrf.m`` script
   (`Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/wrf/matlab/map_wrf.m>`__ or `local
   file <../../models/wrf/matlab/map_wrf.m>`__) can finally plot WRF fields now that the required metadata is part of
   the ``Posterior_Diag.nc``, ``Prior_Diag.nc``, and (not required) ``True_State.nc`` files. It's a small step to
   augment this routine to make publication-quality figures of WRF fields.
-  Regression tests for WRF
   WRF test cases for WRF V2 and V3 for CONUS (Continental or Contiguous United States), a Global WRF case, and a Radar
   test case. The data files are on a web server because they are too large to add to the repository. The README files
   in each directory gives instructions on how to download them.
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/wrf/regression>`__ or `local
   file <../../models/wrf/regression>`__.
-  Other New Model Support
   The ``simple_advection`` and ``MITgcm_ocean`` are fully supported in the Matlab diagnostics.
-  Better execution traces
   Optional detailed execution trace messages from filter by setting the namelist variable ``trace_execution``. See the
   details of the filter namelist
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/filter/filter.html>`__ or
   :doc:`../../assimilation_code/programs/filter/filter` .
-  ``input.nml`` contents saved
   The contents of the ``input.nml`` namelist file are now preserved in the ``True_State.nc``, ``Prior_Diag.nc``, and
   ``Posterior_Diag.nc`` diagnostic files in variable ``inputnml``.
-  Better error checking in obs_sequence creation subroutines to avoid out-of-time-order observations being inserted by
   incorrect programs.
-  Better error checking in ``open_file()``
   Better error checking in the ``utilities_mod`` subroutine ``open_file()``. See documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/modules/utilities/utilities_mod.html#open_file>`__
   or `local file <../../assimilation_code/modules/utilities/utilities_mod.html#open_file>`__.
-  In the DART code tree, individual html pages have links back to the index page, the namelists are moved up to be more
   prominent, and have other minor formatting improvements.
-  The following Matlab observation-space diagnostic routines have been **removed**:

   ======================= ===============================================================================
   fit_ens_mean_time.m     plotted the temporal evolution of the ensemble mean of some quantity.
   fit_ens_spread_time.m   plotted the temporal evolution of the ensemble spread of some quantity.
   fit_mean_spread_time.m  plotted the temporal evolution of the mean and spread of some quantity.
   obs_num_time.m          plotted the temporal evolution of the observation density.
   fit_ens_mean_vertical.m plotted the vertical profile of the ensemble mean of some quantity.
   fit_ens_bias_vertical.m plotted the vertical profile of the bias of the ensemble mean of some quantity.
   obs_num_vertical.m      plotted the vertical profile of the observation density.
   ======================= ===============================================================================

-  The following Matlab observation-space diagnostic routines have been **added**:

   +-----------------------------+---------------------------------------------------------------------------------------+
   | plot_profile.m              | plots the vertical profile of any quantity for any copy with an overlay of the        |
   |                             | observation density and number of observations assimilated. Matlab source             |
   |                             | `Website <https://s                                                                   |
   |                             | vn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/plot_profile.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/plot_profile.m>`__.                          |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | plot_rmse_xxx_profile.m     | plots the vertical profile of the rmse and any quantity for any copy with an overlay  |
   |                             | of the observation density and number of observations assimilated. Matlab source      |
   |                             | `Website <https://svn-dares-                                                          |
   |                             | dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/plot_rmse_xxx_profile.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/plot_rmse_xxx_profile.m>`__.                 |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | plot_bias_xxx_profile.m     | plots the vertical profile of the bias and any quantity for any copy with an overlay  |
   |                             | of the observation density and number of observations assimilated. Matlab source      |
   |                             | `Website <https://svn-dares-                                                          |
   |                             | dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/plot_bias_xxx_profile.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/plot_bias_xxx_profile.m>`__.                 |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | two_experiments_profile.m   | plots the vertical profile of any quantity for any copy for multiple experiments with |
   |                             | an overlay of the observation density and number of observations assimilated in each  |
   |                             | experiment. Matlab source                                                             |
   |                             | `Website <https://svn-dares-da                                                        |
   |                             | rt.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/two_experiments_profile.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/two_experiments_profile.m>`__.               |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | plot_evolution.m            | plots the temporal evolution of any quantity for any copy with an overlay of the      |
   |                             | observation density and number of observations assimilated. Matlab source             |
   |                             | `Website <https://svn                                                                 |
   |                             | -dares-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/plot_evolution.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/plot_evolution.m>`__.                        |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | plot_rmse_xxx_evolution.m   | plots the temporal evolution of the rmse and any quantity for any copy with an        |
   |                             | overlay of the observation density and number of observations assimilated. Matlab     |
   |                             | source                                                                                |
   |                             | `Website <https://svn-dares-da                                                        |
   |                             | rt.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/plot_rmse_xxx_evolution.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/plot_rmse_xxx_evolution.m>`__.               |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | two_experiments_evolution.m | plots the temporal evolution for any quantity for any copy for multiple experiements  |
   |                             | with an overlay of the observation density and number of observations assimilated in  |
   |                             | each experiment. Matlab source                                                        |
   |                             | `Website <https://svn-dares-dart                                                      |
   |                             | .cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/two_experiments_evolution.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/two_experiments_evolution.m>`__.             |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | read_obs_netcdf.m           | reads a netCDF format observation sequence file. Simply need a single copy and a      |
   |                             | single qc - no actual observation required. Matlab source                             |
   |                             | `Website <https://svn-                                                                |
   |                             | dares-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/read_obs_netcdf.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/read_obs_netcdf.m>`__.                       |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | plot_obs_netcdf.m           | reads and plots the locations and values of any copy of the observations in a DART    |
   |                             | netCDF format observation sequence file. Matlab source                                |
   |                             | `Website <https://svn-                                                                |
   |                             | dares-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/plot_obs_netcdf.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/plot_obs_netcdf.m>`__.                       |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | plot_obs_netcdf_diffs.m     | reads and plots the locations and the difference of any two copies of the             |
   |                             | observations in a DART netCDF format observation sequence file. Matlab source         |
   |                             | `Website <https://svn-dares-                                                          |
   |                             | dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/plot_obs_netcdf_diffs.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/plot_obs_netcdf_diffs.m>`__.                 |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | plot_wind_vectors.m         | reads and plots the wind vectors of the observations in a DART netCDF format          |
   |                             | observation sequence file (created by ``obs_seq_to_netcdf``, documentation            |
   |                             | `Website <https://svn-dares-dart.cgd.ucar.edu/DART/re                                 |
   |                             | leases/Kodiak/assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.html>`__ |
   |                             | or :doc:`../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf`)       |
   |                             | Matlab source                                                                         |
   |                             | `Website <https://svn-da                                                              |
   |                             | res-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/plot_wind_vectors.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/plot_wind_vectors.m>`__.                     |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | link_obs.m                  | data brushing tool. Explores many facets of the observations simultaneously. Multiple |
   |                             | plots allow groups of observations to be selected in one view and the corresponding   |
   |                             | observations are indicated in all the other views. Matlab source                      |
   |                             | `Website <https                                                                       |
   |                             | ://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/link_obs.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/link_obs.m>`__.                              |
   +-----------------------------+---------------------------------------------------------------------------------------+
   | plot_rank_histogram.m       | If the individual ensemble member observation values were output from ``filter``      |
   |                             | (selected by namelist option in the filter namelist) into the ``obs_seq.final`` file, |
   |                             | ``obs_diag`` will create rank histogram information and store it in the               |
   |                             | ``obs_diag_output.nc`` file. ``plot_rank_histogram.m`` will then plot it. There are   |
   |                             | instructions on how to view the results with ``ncview`` or with this Matlab script on |
   |                             | the `DART Observation-space                                                           |
   |                             | Diagnos                                                                               |
   |                             | tics <http://www.image.ucar.edu/DAReS/DART/DART_Documentation.php#obs_diagnostics>`__ |
   |                             | web page. Matlab source                                                               |
   |                             | `Website <https://svn-dare                                                            |
   |                             | s-dart.cgd.ucar.edu/DART/releases/Kodiak/diagnostics/matlab/plot_rank_histogram.m>`__ |
   |                             | or `local file <../../diagnostics/matlab/plot_rank_histogram.m>`__.                   |
   +-----------------------------+---------------------------------------------------------------------------------------+

New utilities
-------------

-  obs_seq_to_netcdf
   Any DART observation sequence may be converted to a netCDF format file. All information in the sequence file is
   preserved EXCEPT for any observations with additional user-added metadata, e.g. Radar obs, GPS RO obs for the
   non-local operator. But all core observation data such as location, time, type, QC, observation value and error will
   be converted. This allows for variety of new diagnostics. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.html>`__
   or :doc:`../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf`.
-  obs_seq_coverage
   A step towards determining what locations and quantities are repeatedly observed during a specific time interval.
   This may be used to determine a network of observations that will be used to verify forecasts. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/obs_sequence/assimilation_code/programs/obs_seq_coverage/obs_seq_coverage.html>`__
   or :doc:`../../assimilation_code/programs/obs_seq_coverage/obs_seq_coverage`.
-  obs_selection
   An optional companion routine to ``obs_seq_coverage``. This thins the observation sequence files to contain just the
   desired set of observations to use in the forecast step. This speeds performance by avoiding the cost of evaluating
   observations that will not be used in the verification. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/obs_sequence/assimilation_code/programs/obs_selection/obs_selection.html>`__
   or :doc:`../../assimilation_code/programs/obs_selection/obs_selection`.
-  obs_seq_verify
   is a companion routine to ``obs_seq_coverage``. This creates a netCDF file with variables that should make the
   calculation of skill scores, etc. easier. It creates variables of the form:
   ``METAR_U_10_METER_WIND(analysisT, stations, levels, copy, nmembers, forecast_lead)``
   Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/obs_sequence/assimilation_code/programs/obs_seq_verify/obs_seq_verify.html>`__
   or :doc:`../../assimilation_code/programs/obs_seq_verify/obs_seq_verify`.
-  Select common observation subsets
   A tool that operates on two (will be extended to more) ``obs_seq.final`` files which were output from two different
   runs of filter. Assumes the same ``obs_seq.out`` input file was used in both cases. Outputs two new
   ``obs_seq.final.new`` files containing only the observations which were assimilated in both experiments. It allows
   for a fair comparision with the diagnostic tools. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/programs/obs_common_subset/obs_common_subset.html>`__
   or :doc:`../../assimilation_code/programs/obs_common_subset/obs_common_subset`.
-  Restart File tool
   Generic tool that works on any DART restart file. It is compiled with the corresponding model_mod which tells it how
   large the state vector is. It can alter the timestamps on the data, add or remove model advance times, split a single
   file into 1-per-ensemble or the reverse, and can be used to convert between ASCII and binary formats. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/programs/restart_file_tool/restart_file_tool.html>`__
   or :doc:`../../assimilation_code/programs/restart_file_tool/restart_file_tool`.
-  Advance Time tool
   A generic utility for adding intervals to a Gregorian calendar date and printing out the new date, including handling
   leap year and month and year rollovers. An earlier version of this program was taken from the WRF distribution. This
   version maintains a similar interface but was completely rewritten to use the DART time manager subroutines to do the
   time computations. It reads from the console/standard input to avoid trying to handle command line arguments in a
   compiler-independent manner, and outputs in various formats depending on what is requested via additional flags.
   Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/programs/advance_time/advance_time.html>`__
   or :doc:`../../assimilation_code/programs/advance_time/advance_time`.
-  WRF observation preprocessor tool
   Observation preprocessor which is WRF aware, contributed by Ryan Torn. Will select obs only within the WRF domain,
   will superob, will select only particular obs types based on the namelist. Source is in the
   ``DART/models/wrf/WRF_DART_utilities`` directory.
-  Closest Member tool
   Used in combination with the new option in filter to output the ensemble mean values in a DART restart file format,
   this tool allows you to select the N *closest* members, where there are multiple choices for how that metric is
   computed. There are also ways to select a subset of the state vector by item kind as returned from the
   ``get_state_meta_data()`` routine from the corresponding model interface code in ``model_mod.f90`` (see subroutine
   documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/model_mod.html#get_state_meta_data>`__ or
   `local file <../../models/model_mod.html#get_state_meta_data>`__) and compute the metric based only on those
   values. Tool documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/programs/closest_member_tool/closest_member_tool.html>`__
   or :doc:`../../assimilation_code/programs/closest_member_tool/closest_member_tool`.
-  Fill Inflation restart file tool
   Simple tool that creates an inflation restart file with constant initial inflation and standard deviation values.
   Often the first step of a multi-step assimilation job differs in the namelist only for how the initial inflation
   values are defined. Running this tool creates the equivalent of an IC file for inflation, so the first job step can
   start from a restart file as all subsequent job steps do and allows the use of a single ``input.nml`` file.
   Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/adaptive_inflate/fill_inflation_restart.html>`__
   or :doc:`../../assimilation_code/programs/fill_inflation_restart/fill_inflation_restart`.
-  Replace WRF fields tool
   WRF-specific tool that copies netCDF variables from one file to another. The field must exist in the target file and
   the data will be overwritten by data from the source file. Field names to be copied can be specified directly in the
   namelist or can be listed in a separate file. Missing fields can be ignored or cause the program to stop with a fatal
   error depending on namelist settings. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/models/wrf/WRF_DART_utilities/replace_wrf_fields.html>`__
   or :doc:`../../models/wrf/WRF_DART_utilities/replace_wrf_fields`.
-  model_mod Verification/Check tool
   Tool to help when creating a new model interface file (usually named ``model_mod.f90``). Calls routines to help with
   debugging. Documentation
   `Website <https://svn-dares-dart.cgd.ucar.edu/DART/releases/Kodiak/assimilation_code/programs/model_mod_check/model_mod_check.html%20models/POP/model_mod_check.html>`__
   or :doc:`../../assimilation_code/programs/model_mod_check/model_mod_check`.

Minor items:

-  Most tools which work with observation sequence files now have a namelist option to specify the input files in one of
   two methods: an explicit list of input obs_seq files, or the name of a file which contains the list of obs_seq files.
-  The ``DART/shell_scripts`` directory contains example scripts which loop over multiple days, in formats for various
   shell syntaxes. They are intended as an example for use in advance_model or job scripts, or observation conversion
   programs contributed by users.

Known problems
--------------

We get an internal compiler error when compiling the ``obs_diag`` program on a Linux machine using the gfortran compiler
version 4.1.2. If you get this error try a newer version of the Gnu compiler tools. We have used 4.3 and 4.4
successfully.

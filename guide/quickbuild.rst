.. _DART build system:

DART build system
=================

The DART build system consists of a ``quickbuild.sh`` script for each 
model or observation converter and several build functions in ``DART/build_templates/``

.. code-block:: text

   buildconvfunctions.sh  
   buildfunctions.sh
   buildpreprocess.sh
 

The DART build process is invoked by running ``quickbuild.sh`` and
consists of the following steps:

#. Find the root of the DART git repository you are running in
#. Source the build functions
#. Parse the arguments to quickbuild.sh
#. Clean any existing .o .mod files
#. Run fixsystem to alter compiler dependent source code
#. Compile preprocess
#. Run preprocess to create an obs_def_mod.f90 and obs_kind_mod.f90 specific
   to the observations you are using in DART.
#. For each program

     * collect the source code needed
     * create a Makefile for the executable
     * compile the source code into the executable


The Makefile is created using ``mkmf``, which maps out the dependencies 
between the source files.  For more
information on mkmf please see the `mkmf
documentation <https://github.com/NOAA-GFDL/mkmf>`__.

There is a quickbuild.sh script is each work directory.
To view the usage information for quickbuild.sh,

.. code-block:: text

   ./quickbuild.sh help


``quickbuild.sh`` may be used to build all programs for a particular model or
observation converter, or can be given a single program as an argument to build. 

For example, you may want to build obs_sequence_tool

.. code-block:: bash

   ./quickbuild.sh obs_sequence_tool


In ``quickbuild.sh`` there are arrays containing the list of programs to build.
In general the arrays will contain all programs needed, but you may want to add
other :ref:`DART programs<DART programs>` to your build list. Edit ``quickbuild.sh``
to add your required program to the appropriate array.

For models there are four arrays in quickbuild.sh:

.. code-block:: text
 
  programs=(
  DART programs that can be compiled with mpi go here
  )
   
  serial_programs=(
  DART programs that do not use mpi go here
  )

  model_programs=(
  Model program that can be compiled with mpi go here
  ) 

  model_serial_programs=(
  model programs that do not use mpi go here
  ) 

For observation converters, there is a single array.

.. code-block:: text

     programs=(
     converter programs
     )


For DART developers
--------------------

If you are developing, and want to iterate over changing code and compiling a 
single program, for example filter, you can use ``quickbuild.sh`` in the following 
way:


.. code-block:: bash

   ./quickbuild.sh filter
   *edit the code*
   make
   *edit the code*
   make


Where does quickbuild.sh look for code?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The source code which gets compiled into DART executables is a function of 

* Core DART code
* Location specific code (threed sphere, threed Cartesian, oned, ...)
* Model/converter specific code
* External libraries
* mpi/null mpi utilities
* obs_def_mod.f90 and obs_kind_mod.f90 (created by preprocess)

The core DART code is collected from the ``DART/assimilation_code/modules/`` directory. 
Note the ``DART/assimilation_code/modules/observations`` directory is excluded from 
the search. This directory contains quantity files which are used as input to 
preprocess rather than being compiled directly.

The model directory and the location module to be used are defined in ``quickbuild.sh``.
For example, the Regional Ocean Model (ROMS) uses the threed_sphere location module.

.. code-block ::

   MODEL=ROMS
   LOCATION=threed_sphere

Similarly, for an observation converter, the converter directory and the location 
module are defined in ``quickbuild.sh``

.. code-block ::

   CONVERTER=MADIS
   LOCATION=threed_sphere

The model/observation converter directory will be searched for .f90 files.

*Additionally* any .f90 files in the work directory where you are running 
``quickbuild.sh`` will be added to the list of source files. .f90 files in
the work directory will take precedence over .f90 files with the same name elsewhere. 
 
To take a look how the .f90 files are collected, look at the ``findsrc`` and 
``findconvsrc`` functions in the following files:

.. code-block:: bash

   DART/build_templates/buildfunctions.sh
   DART/build_templates/buildconvfunctions.sh

When adding new code, be sure to obey the following rules to make sure ``quickbuild.sh``
finds your new code and ignores any code you do not want compiled. 

#.  The {name} of the .f90 file must be the program {name}. 
    For example the source code program called ``red_mist`` must 
    be called ``red_mist.f90``

#.  Any .f90 files that you have in your work directory will take precedence over 
    .f90 files with the same name elsewhere. For example if
    you have an ``assim_tools_mod.f90`` in your work directory, this will be 
    compiled rather than the file 
    ``DART/assimilation_code/modules/assimilation/assim_tools_mod.f90``.  
    
    In the example below, the file ``assim_tools_mod.f90`` from the work 
    directory will be used when compiling the lorenz_96 programs.
    
    .. code-block:: text
    
      DART/models/lorenz_96/work/
                                 |-- quickbuild.sh
                                 |-- assim_tools_mod.f90
      
#. If you have .f90 files that you do **not** want to compile into DART, you will
   need to exclude these files using one of these methods:

   * Put the code outside the directories quickbuild.sh searches, for example in a directory 
     ``DART/exclude/``
   * Explicitly exclude the .f90 files with the EXCLUDE variable in ``quickbuild.sh``
   * Rename the .f90 files, e.g. ``solar_flux.f90`` renamed to ``solar_flux.f90.exclude``


#. for core DART programs, use the following directory structure:

   .. code-block:: text

      DART/assimilation_code/programs/{program_name}/
                                                    |-- {program_name}.f90
                                                    |-- {program_name}.rst
                                                    |-- {program_name}.nml
	 
	
   where {program_name}.rst is the documentation for the program and {program_name}.nml
   is a namelist with default values (if applicable to the program).


#. For observation converters, the program must be in the top level of the converter 
   directory:

   .. code-block:: bash

     DART/observations/obs_converters/{converter}/{program_name}.f90


#. For programs specific to a particular model, the program must be in the model directory. 
   For example programs that are specific to Weather and Research Forecasting 
   model (WRF), must be in the ``DART/models/wrf`` directory.
   Model programs may be in subdirectories as shown in the example below. 
   
   .. code-block:: text
   
     DART/models/{model}/
                        |-- {program_one}.f90
                        |-- subdirectory/{program_two}.f90
                               
   
   
   There may be code in the model directory that you do not want compiled into
   the DART executables. For example, the bgrid_solo model directory has all the .f90 code required 
   to build bgrid_solo model (fms_src), which we do *not* want to compile into DART and
   so the bgrid_solo ``quickbuild.sh`` has the following line:
   
   .. code-block:: text
   
      EXCLUDE=fms_src
   
   EXCLUDE is a directory of code to exclude.


#. For code that is outside of the above locations, you can use the variable EXTRA to add source
   files to be compiled.  For example, the ROMS observation converter requires the ROMS 
   model_mod.f90 code, so the ROMS ``quickbuild.sh`` has 

   .. code-block :: text

      EXTRA="$DART/models/ROMS/model_mod.f90"
    
   EXTRA is source code outside the work directory to include in the build. EXTRA can be
   a directory, a list of files, or a single file.



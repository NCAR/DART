DART build system
=================


To view the usage for ``quickbuild.sh``

.. code-block:: bash

  ./quickbuild.sh help


For DART users
--------------

#. compiles preprocess
#. runs preprocess
#. For each program:
     #. collects the source code needed
     #. creates a Makefile for the executable
     #. compiles the source code into the executable

``quickbuild.sh`` may be used to build all programs for a particular model or
observation converter, or can be given a single program as an argument to build. 

For example, you may want to build obs_tool_X

Adding a program:

In ``quickbuild.sh`` there are arrays 

For models,  ``quickbuild.sh``

.. code-block::
      
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

.. code-block:: 

     programs=(
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
----------------------------------------

Essentially, ``quickbuild.sh`` collects all .f90 files from the
``DART/assimilation_code/modules/`` directory. 

Note the ``DART/assimilation_code/modules/observations`` directory is excluded from 
the search. This directory contains quantity files which are used as input to 
preprocess rather than being compiled directly.

Depending on whether ``quickbuild.sh`` is run in a model directory or an 
observation directory additional model/converter directories will be searched.

Additionally any .f90 files in the work directory where you are running 
``quickbuild.sh`` these will be added to the list of source files.
 
  
To take a look how the .f90 files are collected, look at the ``findsrc`` and 
``findconvsrc`` functions in the following files:

.. code-block:: bash

   DART/build_templates/buildfunctions.sh
   DART/build_templates/buildconvfunctions.sh

Important facts:

#.  The {name} of the .f90 file must be the program {name}. 
    For example the source code program called ``red_mist`` must 
    be called ``red_mist.f90``

#. 	Any .f90 files that you have in your work directory will take precedence over 
    other .f90 files. For example if
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
   need to exclude these files.  

   * Put the code outside the directories quickbuild.sh searches, for example in a directory 
     ``DART/exclude/``
   * Explicitly exclude the .f90 files with the EXCLUDE variable in ``quickbuild.sh``
   * Rename the .f90 files, e.g. ``pie.f90`` renamed to ``pie.f90.exclude``


Rules for core DART programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Program directory structure must follow:

.. code-block:: text

   DART/assimilation_code/programs/{program_name}/
                                                 |-- {program_name}.f90
                                                 |-- {program_name}.rst
                                                 |-- {program_name}.nml
	 
	
where {program_name}.rst is the documentation for the program and {program_name}.nml
is a namelist with default values (if applicable to the program).


Rules for observation converter programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a particular converter, the program must be in the top level of the converter 
directory:

.. code-block:: bash

  DART/observations/obs_converters/{converter}/{program_name}.f90


Rules for model programs
~~~~~~~~~~~~~~~~~~~~~~~~~

Programs specific to a particular model must be in the model directory. 
For example programs that are specific to Weather and Research Forecasting 
model (WRF), must be in the ``DART/models/wrf`` directory.
Model programs may be in subdirectories as shown in the example below. 

.. code-block:: text

  DART/models/{model}/
                     |-- {program_one}.f90
                     |-- subdirectory/{program_two}.f90
                                                  
Adding code not in a model directory.

Excluding code in a model directory.









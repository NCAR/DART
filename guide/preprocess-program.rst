How DART supports different types of observations: the preprocess program
=========================================================================

DART's :doc:`preprocess program <../assimilation_code/programs/preprocess/preprocess>`
actually writes the source code that supports
observations. This source code is then used by other modules.

The rationale for ``preprocess``
--------------------------------

Certain types of data require additional metadata in order to be assimilated.
For example, while radiosondes only require the observation location in order
to be assimilated, radar observations need extra metadata to specify the
location of the radar in addition to the location of the observation. GPS
occultations need the locations of the two satellites so the forward operator
can integrate along the raypath. Cosmic ray soil moisture sensors have forward
operators that require site-specific calibration parameters that are not part
of the model and must be included in the observation metadata.

The potential examples are numerous. 

Since each 'observation quantity' may require different amounts of metadata to
be read or written, any routine to read or write an observation sequence
**must** be compiled with support for those particular observations. This is
the rationale for the inclusion of ``preprocess`` in DART. The supported
observations are listed in the ``obs_kind_nml`` namelist in ``input.nml``.

For this reason, we strongly recommend that you use the DART routines to read
and process DART observation sequence files.

.. important::

   You **must** actually run ``preprocess`` before building any executables.
   It is an essential part of DART that enables the same code to interface with
   multiple models and observation types. For example, ``preprocess`` allows
   DART to assimilate synthetic observations for the Lorenz_63 model and real
   radar reflectivities for WRF without needing to specify a set of radar
   operators for the Lorenz_63 model.

``preprocess`` combines multiple ``obs_def`` modules into one
``obs_def_mod.f90`` that is then used by the rest of DART. Additionally, a new
``obs_kind_mod.f90`` is built that will provide support for associating the
specific observation **TYPES** with corresponding (generic) observation
**QUANTITIES**.

The list of ``obs_def`` module source codes is contained in the
``&preprocess_nml`` namelist in ``input.nml``. These modules determine what
observations and operators are supported.

.. warning::
   
   If you want to add another ``obs_def`` module, you **must** rerun
   ``preprocess`` and recompile the rest of your project. ``preprocess``

Example ``preprocess`` namelist
-------------------------------

As an example, if a ``preprocess_nml`` namelist in ``input.nml`` looks like:

.. code-block:: fortran

   &preprocess_nml
       input_obs_kind_mod_file  = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90'
       output_obs_kind_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90'
       quantity_files           = '../../../assimilation_code/modules/observations/atmosphere_quantities_mod.f90',
       input_obs_def_mod_file   = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90'
       obs_type_files           = '../../../observations/forward_operators/obs_def_gps_mod.f90',
                                  '../../../observations/forward_operators/obs_def_QuikSCAT_mod.f90',
                                  '../../../observations/forward_operators/obs_def_GWD_mod.f90',
                                  '../../../observations/forward_operators/obs_def_altimeter_mod.f90',
                                  '../../../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90'
       output_obs_def_mod_file  = '../../../observations/forward_operators/obs_def_mod.f90'
       /

``preprocess`` will combine the following modules:

- ``DEFAULT_obs_def_mod.F90``
- ``obs_def_gps_mod.f90``
- ``obs_def_QuikSCAT_mod.f90``
- ``obs_def_GWD_mod.f90``
- ``obs_def_altimeter_mod.f90``
- and ``obs_def_reanalysis_bufr_mod.f90``
  
into ``obs_def_mod.f90``. This resulting module can be used by the rest of the
project.

Building and running ``preprocess``
-----------------------------------

Since ``preprocess`` is an executable, it must be compiled following the
procedure of all DART executables:

1. The ``DART/build_templates/mkmf.template`` must be correct for your
   environment.
2. The ``preprocess_nml`` namelist in ``input.nml`` must be set properly with
   the modules you want to use.

If those two conditions are met, you can build and run ``preprocess`` using
these commands:

.. code-block::

   $ csh mkmf_preprocess
   $ make
   $ ./preprocess

The first command generates an appropriate ``Makefile`` and the
``input.nml.preprocess_default`` file. The second command results in the
compilation of a series of Fortran90 modules which ultimately produces the
``preprocess`` executable file. The third command actually runs preprocess -
which builds the new ``obs_kind_mod.f90`` and ``obs_def_mod.f90`` source code
files. Once these source code files are created, you can now build the rest of
DART.

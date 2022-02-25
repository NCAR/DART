MODULE DEFAULT_obs_def_mod
==========================

Overview
--------

| DEFAULT_obs_def.F90 is a template used by the program ``preprocess`` to create ``obs_def_mod.f90``.
| To read more detailed instructions on how to add new observation types, see the documentation for
  :doc:`./obs_def_mod`. ``obs_def_*_mod.f90`` files are specified as input to the ``preprocess`` program by namelist,
  and a new ``obs_def_mod.f90`` file is generated which contains all the selected observation types.
| Information from zero or more special obs_def modules, such as ``obs_def_1d_state_mod.f90`` or
  ``obs_def_reanalyis_bufr_mod.f90``, (also documented in this directory) are incorporated into the
  DEFAULT_obs_def_mod.F90 template by ``preprocess``. If no special obs_def files are included in the preprocess
  namelist, a minimal ``obs_def_mod.f90`` is created which can only support identity forward observation operators. Any
  identity observations on the obs_seq.out file will be assimilated, regardless of the obs types specified in
  assimilate_these_obs_types.
| The documentation below describes the special formatting that is included in the ``DEFAULT_obs_def_mod.F90`` in order
  to guide the ``preprocess`` program.

Up to seven sections of code are inserted into ``DEFAULT_obs_def_mod.F90`` from each of the special
``obs_def_*_mod.f90`` files. The insertion point for each section is denoted by a special comment line that must be
included *verbatim* in ``DEFAULT_obs_def_mod.F90``. These special comment lines and their significance are:

#.  ``! DART PREPROCESS MODULE CODE INSERTED HERE``  

   Some special observation definition modules (see for instance ``obs_def_1d_state_mod.f90``) contain code for
   evaluating forward observation operators, reading or writing special information about an observation definition to
   an obs sequence file, or for interactive definition of an observation. The entire module code section is inserted
   here, so the resulting output file will be completely self-contained. Fortran 90 allows multiple modules to be
   defined in a single source file, and subsequent module code can use previously defined modules, so this statement
   must preceed the rest of the other comment lines.

#. ``! DART PREPROCESS USE FOR OBS_QTY_MOD INSERTED HERE``  

   The quantities available to DART are defined by passing quantity files from
   ``DART/assimilation_code/modules/observations`` to ``preprocess``. Unique integer values for each quantity are
   assigned by ``preprocess`` and the use statements for these entries are inserted here.

#. ``! DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE INSERTED HERE``

   Some special observation definition modules (see for instance ``obs_def_1d_state_mod.f90``) contain code for
   evaluating forward observation operators, reading or writing special information about an observation definition to
   an obs sequence file, or for interactive definition of an observation. The use statements for these routines from the
   special observation definition modules are inserted here.

#. ``! DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF INSERTED HERE``

   Special observation definition modules must contain case statement code saying what to do to evaluate a forward
   observation operator for each observation type that they define. This code is inserted here.

#. ``! DART PREPROCESS READ_OBS_DEF INSERTED HERE``

   Special observation definition modules must contain case statement code saying what to do to read any additional
   information required for each observation type that they define from an observation sequence file. This code is
   inserted here.

#. ``! DART PREPROCESS WRITE_OBS_DEF INSERTED HERE``

   Special observation definition modules must contain case statement code saying what to do to write any additional
   information required for each observation type that they define to an observation sequence file. This code is
   inserted here.

#. ``! DART PREPROCESS INTERACTIVE_OBS_DEF INSERTED HERE``

   Special observation definition modules must contain case statement code saying what to do to interactively create any
   additional information required for each observation type that they define. This code is inserted here.

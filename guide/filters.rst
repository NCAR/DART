Filters
=======

The different types of assimilation algorithms (EAKF, ENKF, Kernel filter, Particle filter, etc.) are determined by the
``&assim_tools_nml:filter_kind`` entry, described in :doc:`../../assimilation_code/modules/assimilation/assim_tools_mod`. Despite having
'filter' in the name, they are assimilation algorithms and so are implemented in ``assim_tools_mod.f90``.
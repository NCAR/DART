MODULE DEFAULT_obs_kind_mod
===========================

Overview
--------

DART provides capabilities to assimilate a multitude of different observation
types. Since most DA applications only need to assimilate a subset of the
observation types that DART is capable of assimilating, the observation types 
supported by the programs in your application are defined when you compile
them. You only need to include the observation types you are interested in.

``DEFAULT_obs_kind_mod.F90`` is the input template file which is read by the
:doc:`/assimilation_code/programs/preprocess/preprocess` to create
:doc:`obs_kind_mod`. Information from zero or more special obs_def modules
(such as :doc:`/observations/forward_operators/obs_def_1d_state_mod`) 
and obs_quantities modules 
(such as ``DART/assimilation_code/modules/observations/oned_quantities_mod.f90``) 
are incorporated into the template provided by DEFAULT_obs_def_kind.

If you don't include any specific obs_def files in the preprocessor namelist, 
``preprocess`` will create a minimal ``obs_kind_mod.f90`` file which can only
support identity forward observation operators.

To add a *new* specific observation type, see the
:doc:`../../../observations/forward_operators/obs_def_mod` documentation.

To add a *new* specific observation quantity, see the
:doc:`obs_kind_mod` documentation.

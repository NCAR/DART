.. _adding-your-observations-to-dart:

Adding your observations to DART
================================

First, you should understand that DART already supports a tremendous variety of
observations. To fully support an observation means that the observation can be
converted from its native format to the DART observation sequence format and
that the observation forward operator is already implemented. Keep in mind that
forward operators are not specific to any one model.

:ref:`observations` and :ref:`available_observation_converters`
describe types of observations; their sources and methods for formatting them 
for use by DART. DART can use both real and "synthetic" observations. 
Synthetic observations are used in OSSEs and are useful for model interface development.


Real observations
-----------------

The observation converters are in the *observations/obs_converter* directory,
and are documented in :ref:`available_observation_converters`. 

The forward operators are functionally or logically grouped into Fortran modules
in the *observations/forward_operator* directory. DART employs a 'contractual'
style of programming in that the forward operator requests information from the
model, and if the model cannot provide it, the forward operator may request
different information in an attempt to collect the information needed to apply
the operator. If the model cannot provide any of the required information, the
forward operator fails, the DART QC for that observation is set to the
appropriate value, and the program continues.


Synthetic Observations
-----------------------

Synthetic observations can be created manually using 
:ref:`create_obs_sequence` or extracted from a numerical model using :ref:`pmo`.

When using :ref:`perfect_model_obs<pmo>` the value of the synthetic observation is 
calculated from the model state at a chosen location using the appropriate
forward operator H(x). The observation error variance is specified for later
use and an observation error is added to the value; typically a random draw 
from a distribution sized according to the error variance.

When the observation is used to sample a model state directly, this is known
as an identity observation. The observation error variance and observation 
error are specified in the same way as other synthetic observations, but the 
observation value is a copy of a chosen state variable. That is, the forward 
operator H(x) is the identity matrix (hence the name). The value of the 
observation is taken from the model state at a chosen grid point, which is 
the location of the observation.

Identity observations do not get listed as a type in the header of an 
:ref:`observation sequence file <detailed_structure_obs_seq>`, they are 
denoted in a given observation by a type of -x (negative x) where x is the index in the 
DART state vector that the observation corresponds to. This means that 
identity observations from one model can't (usually) be used in an 
assimilation with another model.

The process for creating synthetic observations, and identity observations is outlined in
:ref:`creating-obs-seq-synthetic` and 
:ref:`synthetic_observations`.

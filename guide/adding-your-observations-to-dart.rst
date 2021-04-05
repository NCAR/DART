Adding your observations to DART
================================

First, you should understand that DART already supports a tremendous variety of
observations. To fully support an observation means that the observation can be
converted from its native format to the DART observation sequence format and
that the observation forward operator is already implemented. Keep in mind that
forward operators are not specific to any one model.

The observation converters are in the *observations/obs_converter* directory and
you should look there for the documentation describing which converters are
available.

The forward operators are functionally or logically grouped into Fortran modules
in the *observations/forward_operator* directory. DART employs a ‘contractual’
style of programming in that the forward operator requests information from the
model, and if the model cannot provide it, the forward operator may request
different information in an attempt to collect the information needed to apply
the operator. If the model cannot provide any of the required information, the
forward operator fails, the DART QC for that observation is set to the
appropriate value, and the program continues.

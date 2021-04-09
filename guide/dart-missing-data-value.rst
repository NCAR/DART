DART missing data value
=======================

If all the prior and posterior mean values are -888888.0 (which is the DART
"missing data" value), those observations were not assimilated.

.. note::
  
   Some observations have precomputed values and the posterior values for these
   will always be -888888.0, no matter if the observation was assimilated or
   not.
   
If it is not already set, edit the ``&filter_nml`` name list in ``input.nml``
to set ``num_output_obs_members`` to be the same as the ensemble size.

This will give you all the forward operator values for all the ensemble
members. You can determine if all ensemble members are failing in the same way,
or if only a few are problematic.

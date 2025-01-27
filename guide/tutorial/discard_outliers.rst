Discarding outlier observations
===============================

In real assimilation applications, it is often useful to discard observations that are too 
far from the model's first guess. 

The quality_control_nml includes outlier_threshold.
Define D as the distance between the prior ensemble mean estimate of an observation and the observation.

If the outlier_threshold is set to x, the observation is not used if D is greater than 
x times the expected separation (see section 5 slide 2).

The default value of outlier_threshold for Lorenz_96 is -1. A negative value means the 
threshold is not applied.

Set the outlier_threshold to 3.
Then rerun filter followed by obs_diag.

Then try rerunning the time series observation space matlab scripts. 
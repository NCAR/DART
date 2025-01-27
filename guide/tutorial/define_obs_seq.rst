DART Observation Sequence Files: Defining Observing Systems
===========================================================

DART uses 'observation sequence' files to define observing systems

The default observation sequence for lorenz_96 has:

	- Observations every 'hour' for 1000 hours,
	- Forty observation stations, randomly located in the periodic 	spatial domain,
	- Observations from each station are taken at each time,
	- The forward operator linearly interpolates from model gridpoints 	to the station location,
	- A random draw from Normal(0, 1) is added to each forward 	operator to simulate observation errors.
    
Five input identity observation sequences that duplicate networks available in DART_LAB guis 
are available in lorenz_96/work:

- obs_seq.in: The default, 40 randomly located observing stations,
- obs_seq_identity.in: Each of the 40 state variables is observed,
- obs_seq_identity_1_40_2.in: Every other state variable is observed,
- obs_seq_identity_1_40_4.in: Every 4th state variable is observed,
- obs_seq_identity_1_20.in: The first 20 state variables are observed.

In all cases, observations are taken every hour for 1000 hours and the observational error 
is randomly selected from Normal(0, 1). 
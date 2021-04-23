Creating an obs_seq file of synthetic observations
==================================================

There are several steps to create an observation sequence file, which follows
directly from the modular nature of the DART programming philosophy.
This procedure may be used to create synthetic observations from *any model*.

1. Decide what observations you want to investigate and edit the
   ``input.nml&obs_kind_nml`` block.
2. Build and run *preprocess* to create code that supports the observations you
   want.
3. Build and run *create_obs_sequence* to define the specifics about the
   observation you want.
4. Build and run *create_fixed_network_sequence* to replicate those specifics
   through time.
5. Build and run *perfect_model_obs* to create an observation consistent with
   the model state and specified error distribution at the requested times and
   locations.

These programs are described in 
:doc:`Programs included in DART <../assimilation_code/programs/readme>`.

Example: generating observations for the Lorenz ’63 model.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While this procedure works with any model, the responses in 'create_obs_sequence'
will vary based on what observations are supported. You should not expect the
responses for observations for L63 can be used to produce radar observations from WRF,
for example. When compiled with support for radar observations, *create_obs_sequence*
**will** prompt you for the required metadata.

1) There are no ‘real’ observations for the Lorenz ’63 model, so the appropriate
namelist settings are:

::

   &obs_kind_nml
       assimilate_these_obs_types = 'RAW_STATE_VARIABLE'  /

   &preprocess_nml
        input_obs_def_mod_file = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90'
       output_obs_def_mod_file = '../../../observations/forward_operators/obs_def_mod.f90'
       input_obs_kind_mod_file = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90'
      output_obs_kind_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90'
                   input_files = '../../../observations/forward_operators/obs_def_1d_state_mod.f90'
     /

2) Run *preprocess* in the normal fashion.

3) *create_obs_sequence* creates an *observation set definition* (typically
named ``set_def.out``), the time-independent part of an observation sequence. It
may help to think of it as trying to define what sorts of observations will be
taken at one ‘reading’ … you walk out to the box and take temperature, humidity,
and wind observations all at the same time and place, for example. You can think
of it as one page in an observer’s notebook, and only contains the *location,
type,* and *observational error characteristics* (normally just the diagonal
observational error variance) for a related set of observations. There are no
actual observation values, nor are there any times associated with the
definition. The program is interactive and queries the user for the information
it needs. Begin by creating a minimal observation set definition in which each
of the 3 state variables of L63 is directly observed with an observational error
variance of 1.0 for each observation. To do this, use the following input
sequence (the text including and after # is a comment and does not need to be
entered):

The following is a screenshot (much of the verbose logging has been left off for
clarity), the user input looks *like this*.

::

      [unixprompt]$ ./create_obs_sequence
       Starting program create_obs_sequence
       Initializing the utilities module.
       Trying to log to unit   10
       Trying to open file dart_log.out

       --------------------------------------
       Starting ... at YYYY MM DD HH MM SS =
                       2017  3 28 10 15 30
       Program create_obs_sequence
       --------------------------------------

       set_nml_output Echo NML values to log file only
       Trying to open namelist log dart_log.nml
       ------------------------------------------------------


       -------------- ASSIMILATE_THESE_OBS_TYPES --------------
       RAW_STATE_VARIABLE
       -------------- EVALUATE_THESE_OBS_TYPES --------------
       ------------------------------------------------------

       ---------- USE_PRECOMPUTED_FO_OBS_TYPES --------------
       ------------------------------------------------------

       Input upper bound on number of observations in sequence
      4
       Input number of copies of data (0 for just a definition)
      0
       Input number of quality control values per field (0 or greater)
      0
       input a -1 if there are no more obs
      0
            Input -1 * state variable index for identity observations
            OR input the name of the observation type from table below:
            OR input the integer index, BUT see documentation...
              1 RAW_STATE_VARIABLE
      -1
       input time in days and seconds
      0 0
       Input error variance for this observation definition
      1.0
       input a -1 if there are no more obs
      0

       { this gets repeated ... until you tell it to stop ... }

       input a -1 if there are no more obs
      -1
       Input filename for sequence (  set_def.out   usually works well)
       set_def.out
       write_obs_seq  opening formatted file set_def.out
       write_obs_seq  closed file set_def.out

Rest assured that if you requested to assimilate more realistic observation
types, you will be queried for appropriate information by *create_obs_sequence*.
Below is a table that explains all of the input you should need to supply for
observations of the L63 model state.

::

   4            # upper bound on num of observations in sequence
   0            # number of copies of data (0 for just a definition)
   0            # number of quality control values per field (0 or greater)
   0            # -1 to exit/end observation definitions

   -1           # observe state variable 1
   0   0        # time -- days, seconds
   1.0          # observational variance
   0            # -1 to exit/end observation definitions

   -2           # observe state variable 2
   0   0        # time -- days, seconds
   1.0          # observational variance
   0            # -1 to exit/end observation definitions

   -3           # observe state variable 3
   0   0        # time -- days, seconds
   1.0          # observational variance
   -1           # -1 to exit/end observation definitions

   set_def.out  # Output file name

4) *create_fixed_network_sequence* takes the observation set definition and
repeats it in time, essentially making multiple pages in our notebook. Again,
the program is interactive and queries the user for information. You should be
able to simply follow the prompts. The table below represents the input needed
for the L63 example:

::

   set_def.out    # Input observation set definition file          
   1              # Regular spaced observation interval in time      
   1000           # 1000 observation times
   0, 43200       # First observation after 12 hours (0 days, 12 * 3600 seconds)
   0, 43200       # Observations every 12 hours
   obs_seq.in     # Output file for observation sequence definition

5) *perfect_model_obs* advances the model from the state defined by the initial
conditions file specified in the ``input.nml`` and ‘applies the forward
operator’ to harvest observations to fill in the observation sequence specified
in ``obs_seq.in``. The observation sequence finally has values for the
observations and is saved in a file generally named *obs_seq.out*.
*perfect_model_obs* is namelist-driven, as opposed to the previous two (whose
input is a lot harder to specify in a namelist). Take a look at (and modify if
you like) the ``input.nml&perfect_model_obs_nml`` section of the namelist.

The End. Not only should you have an observation sequence file (usually
``obs_seq.out``) , you also have a file containing the exact evolution of the
model consistent with those observations - the true state:
``perfect_output.nc``.

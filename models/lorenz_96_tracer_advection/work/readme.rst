.. _quantile tracer: 


Quantile conserving and probit transform tools
==============================================

This file contains instructions for using the lorenz_96_tracer model with DART 
quantile conserving and probit transform filtering tools. These tools are still
being refined, but are working for the examples described. The DART development 
team (dart@ucar.edu) would be happy to hear about your experiences and is
anxious to build scientific collaborations using these new capabilities.


Steps for reproducing basic tests:

Test A: Assimilating observations of state (wind) and tracer concentration using
a rank histogram observation space filter and rank histogram probit transforms for
state variable updates. Example includes adaptive inflation.
The default model configuration has a single tracer source at gridpoint 1 along with
small uniform tracer sinks that lead to areas where the true tracer concentration is
usually 0. This is a particularly tough test for ensemble methods.

#. Build all executables,

   ``./quickbuild.sh nompi`` 
#. Create a set_def.out file using create_obs_sequence:

   ``./create_obs_sequence < create_obs_sequence_input``

#. Create an obs_sequence.in file using ``./create_fixed_network_seq``

   .. code:: text

      ./create_fixed_network_seq
      Select the default input filename <return>,
      Create a regularly repeating sequence by entering "1",
      Enter "1000" for the number of observation times,
      Enter "0 0" for the initial time,
      Enter "0 10800" for the period,
      Select the default output filename, <return>

#. Spin-up a model initial condition by running perfect_model_obs

   ``./perfect_model_obs``

#. Generate a spun-up true time series,

   .. code:: text

      cp perfect_output.nc perfect_input.nc
      Use a text editor to change read_input_state_from_file to .true. in the file input.nml
      Run "./perfect_model_obs" again

#. Run a filter assimilation,

   ``./filter``

#. Examine the output with your favorite tools. Looking at the analysis ensemble 
   for the tracer_concentration variables with indices near the source (location 1)
   and far downstream from the source (location 35) is interesting. Note that the
   source estimation capabilities of the model and filters are not being tested here.


Test B: Using default ensemble adjustment Kalman filters.

The new quantile options are controlled by Fortran code in the module
algorithm_info_mod.f90 in the assimilation_code/modules/assimilation directory.
More information about the control can be found in that module. The tests below 
replace the default version of that module with others that change certain options. 
Doing a diff between these modules shows how the control is being changed for the 
following tests in that module. The tests below 
replace the default version of that module with others that change certain options. 

#. In directory assimilation_code/modules/assimilation, 

   ``cp all_eakf_algorithm_info_mod algorithm_info_mod.f90``

#. Recompile all programs in this directory,

   ``./quickbiuld.sh nompi``

#. Run the filter 
   ``./filter``

Test C: Using default ensemble adjustment Kalman filter for state, but bounded normal rank histogram filter and priors for tracer concentration and source.

#. In directory assimilation_code/modules/assimilation, 

   ``cp state_eakf_tracer_bnrhf_algorithm_info_mod algorithm_info_mod.f90``

#. Recompile all programs in this directory,

   ``./quickbiuld.sh nompi``

#. Run the filter 
   ``./filter``

Test D: Testing bounded above option

Normally tracers are bounded below, but there are other quantities that may be bounded
above. There are distinct numerical challenges in implementing the quantile algorithms
for quantities that are bounded above, so flipping the sign of the tracers is a good
test. 

#. In directory assimilation_code/modules/assimilation,
 
   ``cp neg_algorithm_info_mod algorithm_info_mod.f90``

#. Recompile all programs in this directory,

   ``./quickbiuld.sh nompi``

#. In the file input.nml, change the entry positive_tracer to .false. Also, change the
   entry read_input_state_from_file back to .false. 

#. Repeat steps 3-6 from Test A.



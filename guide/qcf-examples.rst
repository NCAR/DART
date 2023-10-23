.. _quantile tracer:

QCF and Probit Transform Tools: Examples with the Lorenz 96 Tracer Model
========================================================================

This file contains instructions for using the lorenz_96_tracer_advection model with DART 
quantile conserving and probit transform filtering tools. These tools are still
being refined, but are working for the examples described. The DART development 
team (dart@ucar.edu) would be happy to hear about your experiences and is
anxious to build scientific collaborations using these new capabilities.

Make sure that you are on the quantile_methods branch of DART:
``git checkout quantile_methods``

Steps for reproducing basic tests:

Test A: Assimilating observations of state (wind) and tracer concentration using
a rank histogram observation space filter and rank histogram probit transforms for
state variable updates. Example includes adaptive inflation.
The default model configuration has a single tracer source at gridpoint 1 along with
small uniform tracer sinks that lead to areas where the true tracer concentration is
usually 0. This is a particularly tough test for ensemble methods.

#. Download the QCF Table from Google Sheets as a .csv file:
  
   * Visit this link: https://docs.google.com/spreadsheets/d/1ZhKbj0EYKHCgOHvTmJI3k7HI_Ae1NyNKchtekPW0lZs/edit#gid=0
   * Make a copy of the spreadsheet by selecting "File > Make a copy" from the menu bar.
   * Download the spreadsheet as a .csv file by selecting "File > Download > csv" from the menu bar.
   * Google Sheets will append the name of the file with " - Sheet1.csv" when it is downloaded. For example, a spreadsheet named "qcf_table" wil be downloaded as "qcf_table - Sheet1.csv". Rename this file to remove this addition to ensure that there are no spaces in the filename.
   * Copy or move this file to your working directory (/DART/models/lorenz_96_tracer_advection/work).

#. Add the filename of the downloaded .csv file in between the single quotes on the line ``qcf_table_filename = ''`` 
   in the &filter_mod section of /DART/models/lorenz_96_tracer_advection/work/input.nml
   
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

#. Download the QCF Table from Google Sheets as a .csv file:
  
   * Visit this link https://docs.google.com/spreadsheets/d/1e26KuOv_uwrn8y1Ki85FzSeQAc9Pw-nCGk91MpJGVC0/edit#gid=0
   * Make a copy of the spreadsheet by selecting "File > Make a copy" from the menu bar.
   * Download the spreadsheet as a .csv file by selecting "File > Download > csv" from the menu bar.
   * Google Sheets will append the name of the file with " - Sheet1.csv" when it is downloaded. For example, a spreadsheet named "qcf_table" wil be downloaded as "qcf_table - Sheet1.csv". Rename this file to remove this addition to ensure that there are no spaces in the filename.
   * Copy or move this file to your working directory (/DART/models/lorenz_96_tracer_advection/work).

#. Add the filename of the downloaded .csv file in between the single quotes on the line ``qcf_table_filename = ''`` 
   in the &filter_mod section of /DART/models/lorenz_96_tracer_advection/work/input.nml
   
#. Run the filter 
   ``./filter``

Test C: Using default ensemble adjustment Kalman filter for state, but bounded normal rank histogram filter and priors for tracer concentration and source.

#. Download the QCF Table from Google Sheets as a .csv file:
  
   * Visit this link https://docs.google.com/spreadsheets/d/1BEKEnFrw5KI9jf6ewg0POyr98ul5nGjerSVxjqEPDgA/edit#gid=0
   * Make a copy of the spreadsheet by selecting "File > Make a copy" from the menu bar.
   * Download the spreadsheet as a .csv file by selecting "File > Download > csv" from the menu bar.
   * Google Sheets will append the name of the file with " - Sheet1.csv" when it is downloaded. For example, a spreadsheet named "qcf_table" wil be downloaded as "qcf_table - Sheet1.csv". Rename this file to remove this addition to ensure that there are no spaces in the filename.
   * Copy or move this file to your working directory (/DART/models/lorenz_96_tracer_advection/work).

#. Add the filename of the downloaded .csv file in between the single quotes on the line ``qcf_table_filename = ''`` 
   in the &filter_mod section of /DART/models/lorenz_96_tracer_advection/work/input.nml
   
#. Run the filter 
   ``./filter``

Test D: Testing bounded above option

Normally tracers are bounded below, but there are other quantities that may be bounded
above. There are distinct numerical challenges in implementing the quantile algorithms
for quantities that are bounded above, so flipping the sign of the tracers is a good
test. 

#. Download the QCF Table from Google Sheets as a .csv file:
  
   * Visit this link https://docs.google.com/spreadsheets/d/1RHlwyhCpbgcShoQnGW-xp2v-paw1ar-5-EA-uj9CkR8/edit#gid=0
   * Make a copy of the spreadsheet by selecting "File > Make a copy" from the menu bar.
   * Download the spreadsheet as a .csv file by selecting "File > Download > csv" from the menu bar.
   * Google Sheets will append the name of the file with " - Sheet1.csv" when it is downloaded. For example, a spreadsheet named "qcf_table" wil be downloaded as "qcf_table - Sheet1.csv". Rename this file to remove this addition to ensure that there are no spaces in the filename.
   * Copy or move this file to your working directory (/DART/models/lorenz_96_tracer_advection/work).

#. Add the filename of the downloaded .csv file in between the single quotes on the line ``qcf_table_filename = ''`` 
   in the &filter_mod section of /DART/models/lorenz_96_tracer_advection/work/input.nml
   
#. In the file input.nml, change the entry positive_tracer to .false. Also, change the
   entry read_input_state_from_file back to .false. 

#. Repeat steps 5-8 from Test A.

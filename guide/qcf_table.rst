.. _QCF Table:

############################################
Using the QCF Table to Control Input Options
############################################

This file contains instructions for using an input table to set input options with the DART quantile conserving and probit transform filtering tools.
See the following link to learn more about these tools and how to use them:
https://docs.dart.ucar.edu/en/quantile_methods/models/lorenz_96_tracer_advection/work/readme.html

Using this input table allows the user to specify the control options for the Quantile Conserving Filter (QCF), also known as the Quantile Conserving Ensemble Filtering Framework (QCEFF). The observation, state, and inflation variables are all included in this single table.

The new quantile options are read in from the table at runtime and then set in the module algorithm_info_mod.f90 in the DART/assimilation_code/modules/assimilation directory. This module provides routines that give information about details of algorithms for observation error sampling, observation increments, and the transformations for regression and inflation in probit space.

For individual QTYs in DART, the user can specify the options such as the bounds, distribution type, filter kind, etc. for the obs_error_info, probit_dist_info, and obs_inc_info subroutines in algorithm_info_mod.f90

If the user does not use a QCF input table with the DART quantile conserving and probit transform filtering tools, then the default values for these options will be used for all QTYs.

Table Composition
-----------------
Each QTY is specified in its own column, having 28 total control options. 
These control options are divided into 3 main groups, which are the options used for the obs_error_info, probit_dist_info, and obs_inc_info. However, the user is able to specify different values for probit inflation, probit state, and probit extended state, resulting in 5 total groupings for the control options.

The obs_error_info subroutine computes information needed to compute error sample for this observation.
For obs_error_info the input options are the two bounds (lower and upper).

The probit_dist_info subroutine computes the details of the probit transform.
From probit_dist_info, the values needed are the bounds and the distribution type. These can be different for all three cases (inflation, state, and extended_state).

The obs_inc_info subrotuine sets the details of how to assimilate this observation.
From obs_inc_info, the bounds, plus the filter_kind, rectangular_quadrature, gaussian_likelihood_tails, sort_obs_inc, and spread_restoration are needed. However, rectangular_quadrature and gaussian_likelihood_tails are only applicable with RHF.

Full list of options:
Obs_error_info: bounded_below, bounded_above, lower_bound, upper_bound [4 columns]
Probit_dist_info: dist_type, bounded_below, bounded_above, lower_bound, upper_bound (x3 for inflation, state, and observation (extended state) priors) [15 columns]
Obs_inc_info: filter_kind, rectangular_quadrature, gaussian_likelihood_tails, sort_obs_inc, spread_restoration, bounded_below, bounded_above, lower_bound, upper_bound [9 columns]

Customizing the Table
---------------------
The table can be customized by either editing a YAML file (which is then converted to a tabular data file in .txt format by a Python script) or a Google Sheet spreadsheet (which is then downloaded in .csv format). The specifics of how to manually edit both formats will be detailed in the following sections.

Regardless of which of these formats are used, the table consists of two headers. The first states the version # of the table being used; the most recent version of the table needs to be used to ensure compatibilty with DART. The current version # is 1. The second header lists the full set of input options, or all 28 column names in other words.

Generally, the user will add and fill in one row for each bounded QTY. If a QTY is not listed in the table, the default values will be used for all 28 options. Therefore, the user will only need to add rows for QTYs that use non-default values for any of the input options. **

The majority of the input options are read in as logicals, and will need to be written in the format of either 'F' or '.false.' These include bounded_below, bounded_above, rectangular_quadrature, gaussian_likelihood_tails, sort_obs_inc, and spread_restoration.

The actual numerical values of the bounds are read in as real_r8 types. These can be specified as reals or integers in the table. 

dist_type and filter_kind are read in as strings, which the possible values for are listed below:

dist_type:
NORMAL_DISTRIBUTION
BOUNDED_NORMAL_RH_DISTRIBUTION
GAMMA_DISTRIBUTION
BETA_DISTRIBUTION
LOG_NORMAL_DISTRIBUTION
UNIFORM_DISTRIBUTION
PARTICLE_FILTER_DISTRIBUTION

filter_kind:
EAKF
ENKF
UNBOUNDED_RHF
GAMMA_FILTER
BOUNDED_NORMAL_RHF

The default values for each of the options are listed below:
bounded_below = .false.
bounded_above = .false.
lower_bound   = -888888
upper_bound   = -888888
dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
filter_kind = BOUNDED_NORMAL_RHF
rectangular_quadrature = .false.
gaussian_likelihood_tails = .false.
sort_obs_inc = .false.
spread_restoration = .false.

Note that bounds set to -888888 are missing_r8 values.

YAML File Usage
---------------
This section will detail how to customize the qcf_table_template.yaml file and then utilize the yaml_to_table.py Python script to convert the YAML dictionary into a table in .txt format.

First, the user needs to access YAML template file, located in DART/assimilation/programs/qcf_table/
This template file is then to be copied into another file. You can name this anything, but the standard name is 'qcf_table.yaml'.

.. code::
  cp qcf_table_template.yaml qcf_table.yaml

The YAML file needs to match the formatting in qcf_table_template.yaml, which is as follows:

::

   QCF table version: 1
   QTY_TEMPLATE:
     obs_error_info:
       bounded_below
       bounded_above
       lower_bound
       upper_bound
     probit_inflation:
       dist_type
       bounded_below
       bounded_above
       lower_bound
       upper_bound
     probit_state:
       dist_type
       bounded_below
       bounded_above
       lower_bound
       upper_bound
     probit_extended_state:
       dist_type
       bounded_below
       bounded_above
       lower_bound
       upper_bound
     obs_inc_info:
       filter_kind
       rectangular_quadrature
       gaussian_likelihood_tails
       sort_obs_inc
       spread_restoration
       bounded_below
       bounded_above
       lower_bound
       upper_bound
   QTY_STATE_VARIABLE:
     obs_error_info:
       bounded_below: .false.
       bounded_above: .false.
       lower_bound: -888888.0
       upper_bound: -888888.0
     probit_inflation:
       dist_type: BOUNDED_NORMAL_RH_DISTRIBUTION
       bounded_below: .false.
       bounded_above: .false.
       lower_bound: -888888.0
       upper_bound: -888888.0
     probit_state:
       dist_type: BOUNDED_NORMAL_RH_DISTRIBUTION
       bounded_below: .false.
       bounded_above: .false.
       lower_bound: -888888.0
       upper_bound: -888888.0
     probit_extended_state:
       dist_type: BOUNDED_NORMAL_RH_DISTRIBUTION
       bounded_below: .false.
       bounded_above: .false.
       lower_bound: -888888.0
       upper_bound: -888888.0
     obs_inc_info:
       filter_kind: BOUNDED_NORMAL_RHF
       rectangular_quadrature: .false.
       gaussian_likelihood_tails: .false.
       sort_obs_inc: .false.
       spread_restoration: .false.
       bounded_below: .false.
       bounded_above: .false.
       lower_bound: -888888.0
       upper_bound: -888888.0

To customize the YAML dictionary file, the user should change the name 'QTY_STATE_VARIABLE' to the name of the first QTY to be specified with non-default values. Edit the values for the vairables wanting to be changed, and leave the rest of the variables set to the default values.

To add additional QTYs after this, simply copy the lines pertaining to first QTY, change the name of the QTY, and set the variables accordingly.

To remove a QTY from the YAML dictionary, simply remove the lines it consists of.

The user will then take their customized YAML file and pass it as input into a Python script. This will convert it into a text file contaning the table data. 

This script is located in DART/assimilation/programs/qcf_table/

To use the Python script on Derecho or Cheyenne, the user must first load the correct modules

::

   module load conda
   conda activate npl

Then run the python script.

::

   python3 yaml_to_table.py

The user will be prompted to enter the name of the input YAML file and the name for the output text file name.
A table will be produced at the specified output filename.

Copy or move this file to your working directory.

Google Sheets Usage
-------------------
This section will detail how to customize the Google Sheets spreadsheet and then download the spreadsheet into a table in .csv format.

Folow this link https://docs.google.com/spreadsheets/d/1SI4wHBXatLAAMfiMx3mUUC7x0fqz4lniKuM4_i5j6bM/edit#gid=0 to access the template spreadsheet.

The QTYs listed in the template file (QTY_STATE_VARIABLE, QTY_TRACER_SOURCE) correspond to the lorenz_6_tracer_advection model and have the default values set for all variables. Make sure to remove these QTYs if you are not running an analagous model. **

Make a copy of the table by selecting 'File > Make a copy' from the menu bar.

To customize the spreadsheet, click on the cell you want to edit and change the value of that cell.
To add a new QTY to the spreadsheet, simply copy the row of a listed QTY, change the QTY name, and edit the cells individually to set the control options.
To remove a QTY from the spreadsheet, select the row corresponding to that QTY. Then right click and choose "Delete Row"

Ensure that there are no empty rows in between the QTYs listed in the spreadsheet.

Download the spreadsheet as a .csv file by selecting 'File > Download > csv' from the menu bar.

Google Sheets will append the name of the file with " - Sheet1.csv". For example a spreadsheet named "qcf_table" wil be downloaded as "qcf_table - Sheet1.csv" 
Rename this file to remove this append to ensure that there are no spaces in the filename.

Copy or move this file to your working directory.

Using the table in DART
-----------------------
Navigate to your working directory.

Edit your namelist file (input.nml)
Add the item "qcf_table_filename = 'your_filename' to the &filter_nml section, replacing your_filename with the actual name of the file you want to use.
Remember that the default values will be used for all QTYs if no filename is listed here.

Build and run filter normally.

The data read from the QCF table used is written to the output file dart_log.out

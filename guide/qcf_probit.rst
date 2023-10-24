.. _QCF:

########################################################
Quantile Conserving and Probit Transform Filtering Tools
########################################################

This file contains instructions for using the DART Quantile Conserving Filters (QCF), also known as the Quantile Conserving Ensemble Filtering Framework (QCEFF), and probit transform filtering tools.

The DART development team (dart@ucar.edu) would be happy to hear about your experiences and is anxious to build scientific collaborations using these new capabilities.

The user can include an input table allows the user to specify the control options for these tools. The observation, state, and inflation variables are all included in this single table.

The new quantile options are read in from the table at runtime and then set in the module algorithm_info_mod.f90 in the DART/assimilation_code/modules/assimilation directory. This module provides routines that give information about details of algorithms for observation error sampling, observation increments, and the transformations for regression and inflation in probit space.

For individual QTYs in DART, the user can specify the options such as the bounds, distribution type, filter kind, etc. for the obs_error_info, probit_dist_info, and obs_inc_info subroutines in algorithm_info_mod.f90

If the user does not use a QCF input table with the DART quantile conserving and probit transform filtering tools, then the default values for these options will be used for all QTYs.

Table Composition
-----------------
The table consists of two headers. The first states the version # of the table being used; the most recent version of the table needs to be used to ensure compatibilty with DART. The current version # is 1. The second header lists the full set of input options, or all 24 column names in other words.

Each QTY is specified in its own column, having 24 total control options. 
These control options are divided into 3 main groups, which are the options used for the obs_error_info, probit_dist_info, and obs_inc_info. However, the user is able to specify different values for probit inflation, probit state, and probit extended state, resulting in 5 total groupings for the control options.

The obs_error_info subroutine computes information needed to compute error sample for this observation.
For obs_error_info the input options are the two bounds (lower and upper).

The probit_dist_info subroutine computes the details of the probit transform.
From probit_dist_info, the values needed are the bounds and the distribution type. These can be different for all three cases (inflation, state, and extended_state).

The obs_inc_info subrotuine sets the details of how to assimilate this observation.
From obs_inc_info, the values needed are the bounds and the filter_kind.

Full list of options:
Obs_error_info: bounded_below, bounded_above, lower_bound, upper_bound [4 columns]
Probit_dist_info: dist_type, bounded_below, bounded_above, lower_bound, upper_bound (x3 for inflation, state, and observation (extended state) priors) [15 columns]
Obs_inc_info: filter_kind, bounded_below, bounded_above, lower_bound, upper_bound [5 columns]

Customizing the Table
---------------------
The table can be customized by editing a Google Sheet spreadsheet (which is then downloaded in .csv format). Folow this `link <https://docs.google.com/spreadsheets/d/1SI4wHBXatLAAMfiMx3mUUC7x0fqz4lniKuM4_i5j6bM/edit#gid=0>`_ to access the template spreadsheet.

The user will add and fill in one row for each bounded QTY they want to specify. If a QTY is not listed in the table, the default values will be used for all 25 options. Therefore, the user will only need to add rows for QTYs that use non-default values for any of the input options.

The default values for each of the options are listed below:
bounded_below = .false.
bounded_above = .false.
lower_bound   = -888888
upper_bound   = -888888
dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
filter_kind = BOUNDED_NORMAL_RHF

Note that bounds set to -888888 are missing_r8 values.

The following input options are read in as logicals, and will need to be written in the format of either 'F' or '.false.' These include bounded_below, bounded_above, and spread_restoration.

The actual numerical values of the bounds are read in as real_r8 types. These can be specified as reals or integers in the table. 

dist_type and filter_kind are read in as strings. The possible values for these variables are listed below:

dist_type:
NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, GAMMA_DISTRIBUTION, BETA_DISTRIBUTION, LOG_NORMAL_DISTRIBUTION, UNIFORM_DISTRIBUTION, PARTICLE_FILTER_DISTRIBUTION

filter_kind:
EAKF, ENKF, UNBOUNDED_RHF, GAMMA_FILTER, BOUNDED_NORMAL_RHF

Make a copy of the table by selecting 'File > Make a copy' from the menu bar.

To customize the spreadsheet, click on the cell you want to edit and change the value of that cell.
To add a new QTY to the spreadsheet, copy row 3 of the table into the next available row, change ``QTY_NAME`` to the name of the QTY to specify, and edit the cells individually to set the control options.
To remove a QTY from the spreadsheet, select the row number corresponding to that QTY. Then right click and choose "Delete Row"
Make sure to remove the row for ``QTY_NAME`` when you have finished adding all of the specified QTYs to the table.

Ensure that there are no empty rows in between the QTYs listed in the spreadsheet.

Download the spreadsheet as a .csv file by selecting 'File > Download > csv' from the menu bar.

Google Sheets will append the name of the file with " - Sheet1.csv" when it is downloaded. For example, a spreadsheet named "qcf_table" wil be downloaded as "qcf_table - Sheet1.csv" 
Rename this file to remove this addition to ensure that there are no spaces in the filename.

Copy or move this file to your working directory (/DART/models/model_name/work).

Using the table in DART
-----------------------
Navigate to your working directory (/DART/models/model_name/work).

Switch to the quantile_methods branch of DART:
``git checkout quantile_methods``

Edit your namelist file (input.nml):
Add the name of the QCF table file in between the quotes of ``qcf_table_filename = ''`` in the &filter_nml section.
Remember that the default values will be used for all QTYs if no filename is listed here.

Build and run filter normally.

The data that is read from in the QCF table is written to the output file dart_log.out

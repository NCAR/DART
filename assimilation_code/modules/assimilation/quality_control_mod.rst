MODULE quality_control_mod
==========================

Overview
--------

Routines in this module deal with two different types of quality control (QC) related functions. The first is to support
interpretation of the *incoming* data quality, to reject observations at assimilation time which are marked as poor
quality. The second is to document how DART disposed of each observation; whether it was successfully assimilated or
rejected, and if rejected, for which reason.

Usage
-----

Incoming data quality control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DART currently supports a single incoming quality control scheme compatible with NCEP usage. Lower values are considered
better and higher values are considered poorer. A single namelist item, ``input_qc_threshold`` sets the boundary between
accepted and rejected observations. Values *larger* than this value are rejected; values equal to or lower are accepted.
Note that observations could be subsequently rejected for other reasons, including failing the outlier threshold test or
all observations of this type being excluded by namelist control. See the
`obs_kind_mod <../observations/obs_kind_mod.html#Namelist>`__ namelist documentation for more details on how to enable
or disable assimilation by observation type at runtime.

The incoming quality control value is set when an observation sequence file is created. If the data provider user a
different scheme the values must be translated into NCEP-consistent values. Generally we use the value 3 for most runs.

Observations can also be rejected by the assimilation if the observation value is too far from the mean of the ensemble
of expected values (the forward operator results). This is controlled by the ``outlier_threshold`` namelist item.

Specifically, the outlier test computes the difference between the observation value and the prior ensemble mean. It
then computes a standard deviation by taking the square root of the sum of the observation error variance and the prior
ensemble variance for the observation. If the difference between the ensemble mean and the observation value is more
than the specified number of standard deviations then the observation is not used. This can be an effective way to
discard clearly erroneous observation values. A commonly used value is 3. To assimilate all possible observations, 
a value of -1 can be used, but may result in 'chasing bad observations' and _prevents_ the calculation of the number 
of observations that are grossly inconsistent with the ensemble; a useful indicator of *filter divergence*.

There is an option to add code to this module to specialize the outlier threshold routine. For example, it is possible
to allow all observations of one type to be assimilated regardless of the outlier value, and enforce the outlier
threshold only on other types of observations. To enable this capability requires two actions: setting the
``enable_special_outlier_code`` namelist to ``.TRUE.``, and adding your custom code to the ``failed_outlier()``
subroutine in this module.

DART outgoing quality control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As DART assimilates each observation it adds a *DART Quality Control* value to the output observation sequence
(frequently written to a file named ``obs_seq.final)``. This flag indicates how the observation was used during the
assimilation. The flag is a numeric value with the following meanings:

== ====================================================================================================================
0: Observation was assimilated successfully
1: Observation was evaluated only so not used in the assimilation
2: The observation was used but one or more of the posterior forward observation operators failed
3: The observation was evaluated only so not used AND one or more of the posterior forward observation operators failed
4: One or more prior forward observation operators failed so the observation was not used
5: The observation was not used because it was not selected in the namelist to be assimilated or evaluated
6: The incoming quality control value was larger than the threshold so the observation was not used
7: Outlier threshold test failed (as described above)
8: The location conversion to the vertical localization unit failed so the observation was not used
== ====================================================================================================================

| 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &quality_control_nml
      input_qc_threshold          = 3
      outlier_threshold           = -1
      enable_special_outlier_code = .false.
     /

| 

Items in this namelist control whether an observation is assimilated or not.

.. container::

   +-----------------------------+----------+----------------------------------------------------------------------------+
   | Item                        | Type     | Description                                                                |
   +=============================+==========+============================================================================+
   | input_qc_threshold          | real(r8) | Numeric value indicating whether this observation is considered "good      |
   |                             |          | quality" and should be assimilated, or whether it is suspect because of    |
   |                             |          | previous quality control processes. This value would have been set when    |
   |                             |          | the observation was created and added to the observation sequence file.    |
   |                             |          | Observations with an incoming QC value larger than this threshold are      |
   |                             |          | rejected and not assimilated.                                              |
   +-----------------------------+----------+----------------------------------------------------------------------------+
   | outlier threshold           | real(r8) | This numeric value defines the maximum number of standard deviations an    |
   |                             |          | observation value can be away from the ensemble forward operator mean and  |
   |                             |          | still be assimilated. Setting it to the value -1 disables this check.      |
   +-----------------------------+----------+----------------------------------------------------------------------------+
   | enable_special_outlier_code | logical  | Setting this value to .TRUE. will call a subroutine ``failed_outlier()``   |
   |                             |          | instead of using the default code. The user can then customize the tests   |
   |                             |          | in this subroutine, for example to accept all observations of a particular |
   |                             |          | type, or use different numerical thresholds for different observation      |
   |                             |          | types or locations.                                                        |
   +-----------------------------+----------+----------------------------------------------------------------------------+

| 

Discussion
----------

Small ensemble spread
^^^^^^^^^^^^^^^^^^^^^

If an ensemble is spun up from a single state the ensemble spread may be very small to begin and many observations may
be rejected by the ``outlier_threshold``. But as the ensemble spread increases the assimilation should be able to
assimilate more and more observations as the model trajectory becomes consistent with those observations.

Other modules used
------------------

::

   types_mod
   utilities_mod
   random_seq_mod

Public interfaces
-----------------

=================================== =======================
``use quality_control_mod, only :`` initialize_qc
\                                   input_qc_ok
\                                   get_dart_qc
\                                   check_outlier_threshold
\                                   good_dart_qc
\                                   set_input_qc
\                                   dart_flags
=================================== =======================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *call check_outlier_threshold(obs_prior_mean, obs_prior_var, obs_val, obs_err_var, & obs_seq, this_obs_key, dart_qc)*
   ::

      real(r8),                intent(in)    :: obs_prior_mean !>  prior observation mean
      real(r8),                intent(in)    :: obs_prior_var  !>  prior observation variance
      real(r8),                intent(in)    :: obs_val        !>  observation value
      real(r8),                intent(in)    :: obs_err_var    !>  observation error variance
      type(obs_sequence_type), intent(in)    :: obs_seq        !>  observation sequence
      integer,                 intent(in)    :: this_obs_key   !>  index for this observation
      integer,                 intent(inout) :: dart_qc        !>  possibly modified DART QC

.. container:: indent1

   Computes whether this observation failed the outlier threshold test and if so, updates the DART QC.

| 

.. container:: routine

   *var = input_qc_ok(input_qc, qc_to_use)*
   ::

      real(r8), intent(in)  :: input_qc    !> incoming QC data value
      integer,  intent(out) :: qc_to_use   !> resulting DART QC
      logical               :: input_qc_ok !> true if input_qc is good

.. container:: indent1

   Returns true if the input qc indicates this observation is good to use.

| 

.. container:: routine

   ::

      ! Dart quality control variables
      integer, parameter :: DARTQC_ASSIM_GOOD_FOP        = 0
      integer, parameter :: DARTQC_EVAL_GOOD_FOP         = 1
      integer, parameter :: DARTQC_ASSIM_FAILED_POST_FOP = 2
      integer, parameter :: DARTQC_EVAL_FAILED_POST_FOP  = 3
      integer, parameter :: DARTQC_FAILED_FOP            = 4
      integer, parameter :: DARTQC_NOT_IN_NAMELIST       = 5
      integer, parameter :: DARTQC_BAD_INCOMING_QC       = 6
      integer, parameter :: DARTQC_FAILED_OUTLIER_TEST   = 7
      integer, parameter :: DARTQC_FAILED_VERT_CONVERT   = 8
      !!integer, parameter :: DARTQC_OUTSIDE_DOMAIN        = 9  ! we have no way (yet) for the model_mod to signal this

.. container:: indent1

   These are public constants for use in other parts of the DART code.

| 

Files
-----

========= ========================================
filename  purpose
========= ========================================
input.nml to read the quality_control_mod namelist
========= ========================================

References
----------

#. none

Error codes and conditions
--------------------------

.. container:: errors

   ============ ============= ======================
   Routine      Message       Comment
   ============ ============= ======================
   routine name output string description or comment
   ============ ============= ======================

Future plans
------------

Should support different incoming data QC schemes.

It would be nice to have a different DART QC flag for observations which fail the forward operator because they are
simply outside the model domain. The diagnosic routines may indicate a large number of failed forward operators which
make it confusing to identify observations where the forward operator should have been computed and can skew the
statistics. Unfortunately, this requires adding an additional requirement on the model-dependent *model_mod.f90* code in
the ``model_interpolate()`` routine. The current interface defines a return status code of 0 as success, any positive
value as failure, and negative numbers are reserved for other uses. To identify obs outside the domain would require
reserving another value that the interpolate routine could return.

At this time the best suggestion is to cull out-of-domain obs from the input observation sequence file by a
preprocessing program before assimilation.

Private components
------------------

N/A

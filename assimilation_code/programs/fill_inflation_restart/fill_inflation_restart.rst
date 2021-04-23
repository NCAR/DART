PROGRAM ``fill_inflation_restart``
==================================

Overview
--------

Utility program to create inflation restart files with constant values.

These files can be used as input for the first step of a multi-step assimilation when adaptive inflation is being used.
This allows the namelist items ``inf_initial_from_restart`` and ``inf_sd_initial_from_restart`` in the ``&filter_nml``
namelist to be ``.TRUE.`` for all steps of the assimilation including the very first one. (These items control whether
inflation values are read from an input file or read from constants in the namelist.)

Adaptive inflation restart files are written at the end of a ``filter`` run and are needed as input for the next
timestep. This program creates files that can be used for the initial run of filter when no inflation restart files have
been created by filter but are required to be read as input.

This program reads the inflation values to use from the ``&fill_inflation_restart_nml`` namelist for setting the prior
inflation mean and standard deviation, and/or the posterior inflation mean and standard deviation. It does not use the
inflation values in the ``&filter`` namelist.

This program uses the information from the model_mod code to determine the number of items in the state vector. It must
be compiled with the right model's model_mod, and if the items in the state vector are selectable by namelist options,
the namelist when running this program must match exactly the namelist used during the assimilation run.

It creates files with names consistent with the input names expected by filter:

::

   input_priorinf_mean.nc
   input_priorinf_sd.nc
   input_postinf_mean.nc
   input_postinf_sd.nc

An older (and deprecated) alternative to running ``fill_inflation_restart`` is to create inflation netcdf files by using
one of the NCO utilities like "``ncap2``" on a copy of another restart file to set the initial inflation mean, and
another for the initial inflation standard deviation. Inflation mean and sd values look exactly like restart values,
arranged by variable type like T, U, V, etc.

Depending on your version of the NCO utilities, you can use ``ncap2`` to set the T,U and V inf values using one of two
syntaxes:

.. container:: unix

   ::

        ncap2 -s 'T=1.0;U=1.0;V=1.0' wrfinput_d01 input_priorinf_mean.nc
        ncap2 -s 'T=0.6;U=0.6;V=0.6' wrfinput_d01 input_priorinf_sd.nc
        -or-
        ncap2 -s 'T(:,:,:)=1.0;U(:,:,:)=1.0;V(:,:,:)=1.0' wrfinput_d01 input_priorinf_mean.nc
        ncap2 -s 'T(:,:,:)=0.6;U(:,:,:)=0.6;V(:,:,:)=0.6' wrfinput_d01 input_priorinf_sd.nc

Some versions of the NCO utilities change the full 3D arrays into a single scalar. If that's your result (check your
output with ``ncdump -h``) use the alternate syntax or a more recent version of the NCO tools.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &fill_inflation_restart_nml

      write_prior_inf    = .FALSE.
      prior_inf_mean     = -88888.8888
      prior_inf_sd       = -88888.8888

      write_post_inf     = .FALSE.
      post_inf_mean      = -88888.8888
      post_inf_sd        = -88888.8888

      single_file        = .FALSE.
      input_state_files  = ''
      verbose            = .FALSE.
   /

The namelist controls which files are created and what values are written to the restart files.

.. container::

   +-------------------+--------------+---------------------------------------------------------------------------------+
   | Item              | Type         | Description                                                                     |
   +===================+==============+=================================================================================+
   | write_prior_inf   | logical      | Setting this to .TRUE. writes both the prior inflation mean and standard        |
   |                   |              | deviation files: ``input_priorinf_mean.nc``, ``input_priorinf_sd.nc``.          |
   +-------------------+--------------+---------------------------------------------------------------------------------+
   | prior_inf_mean    | real(r8)     | Prior inflation mean value.                                                     |
   +-------------------+--------------+---------------------------------------------------------------------------------+
   | prior_inf_sd      | real(r8)     | Prior inflation standard deviation value.                                       |
   +-------------------+--------------+---------------------------------------------------------------------------------+
   | write_post_inf    | logical      | Setting this to .TRUE. writes both the posterior inflation mean and standard    |
   |                   |              | deviation files ``input_postinf_mean.nc``, ``input_postinf_sd.nc``.             |
   +-------------------+--------------+---------------------------------------------------------------------------------+
   | post_inf_mean     | real(r8)     | Posterior inflation mean value.                                                 |
   +-------------------+--------------+---------------------------------------------------------------------------------+
   | post_inf_sd       | real(r8)     | Posterior inflation standard deviation value.                                   |
   +-------------------+--------------+---------------------------------------------------------------------------------+
   | single_file       | logical      | Currently not supported, but would be used in the case where you have a single  |
   |                   |              | restart file that contains all of the ensemble members. Must be .false.         |
   +-------------------+--------------+---------------------------------------------------------------------------------+
   | input_state_files | character(:) | List one per domain, to be used as a template for the output inflation files.   |
   +-------------------+--------------+---------------------------------------------------------------------------------+
   | verbose           | logical      | Setting this to .TRUE. gives more output, and is generally used for debugging   |
   +-------------------+--------------+---------------------------------------------------------------------------------+

| 

Here is an example of a typical namelist for ``fill_inflation_restart`` :

::

   &fill_inflation_restart_nml

      write_prior_inf    = .TRUE.
      prior_inf_mean     = 1.01
      prior_inf_sd       = 0.6

      write_post_inf     = .FALSE.
      post_inf_mean      = 1.0
      post_inf_sd        = 0.6

      single_file        = .FALSE.
      input_state_files  = ''
      verbose            = .FALSE.
   /

Files
-----

Creates:

::

   input_priorinf_mean.nc
   input_priorinf_sd.nc
   input_postinf_mean.nc
   input_postinf_sd.nc

based on the template file from the specific model this code is compiled for.

References
----------

-  none

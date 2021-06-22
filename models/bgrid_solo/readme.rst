bgrid_solo
==========

Overview
--------

DART interface module for the dynamical core of the GFDL AM2 Bgrid model. This
model is subroutine callable from DART and can be run in a similar fashion to
low-order models that produce diagnostic output files with multiple assimilation
times per file.

The Bgrid model was originally configured as a comprehensive atmospheric model
as described in Anderson et al. (2004). [1]_

All of that code remains in the directories under the
``DART/models/bgrid_solo`` directory, however, much of the capability has
been disabled by code modification. What is left is a dry dynamical core for a
model with no diurnal cycle at equinox with forcing described in Held and Suarez
(1994). [2]_

The default settings are for a model with a 60x30 horizontal grid and 5 vertical
levels. This is close to the smallest version that has somewhat realistic
baroclinic instability resulting in mid-latitude 'storm tracks'. The model
resolution can be changed with the entries in the ``bgrid_cold_start_nml``
namelist described in the `Namelist`_ section. It may be necessary to change the
model time step to maintain stability for larger model grids. The model state
variables are the gridded surface pressure, temperature, and u and v wind
components.

The ``bgrid_solo`` directory has a ``work/workshop_setup.csh`` script that compiles 
and runs an example. This example is intended to demonstrate that the same
process used for a low-order model may be used for a much more 
complex model and generates output for state-space or observation-space diagnostics. 

Some examples of ways in which this model can be configured and modified to test
DART assimilation capabilities are documented in Anderson et al. (2005). [3]_

Several programs that generate interesting observation sequences are available
in the ``DART/models/bgrid_solo`` directory. These programs take
interactive user input and create a text file that can be piped into program
``create_obs_sequence`` to create obs_sequence files. These can serve as
examples for users who are interested in designing their own custom obs_sequence
files.

Program ``column_rand`` creates an obs_sequence with randomly located columns of
observations (essentially synthetic radiosondes) that observe surface pressure
along with temperature and wind components at all model levels.

Program ``id_set_def_stdin`` generates an obs_sequence file that observes every
state variable with error variance of 10000 for surface pressure and 1.0 for
temperature and wind components.

Program ``ps_id_stdin`` generates an obs_sequence that observes every surface
pressure variable for the default model size (30x60) with an error variance of
100.

Program ``ps_rand_local`` generates a set of randomly located surface pressure
observations with an interactively specified error variance. It also allows the
observations to be confined to a rectangular subdomain.

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.
 
.. code-block:: fortran

   &model_nml 
      current_time =  0, 0, 0, 0
      override = .false.,
      dt_atmos = 3600,
      days     = 10,
      hours    = 0,
      minutes  = 0,
      seconds  = 0,
      noise_sd = 0.0,
      dt_bias  = -1,
      state_variables = 'ps', 'QTY_SURFACE_PRESSURE',
                        't',  'QTY_TEMPERATURE',
                        'u',  'QTY_U_WIND_COMPONENT',
                        'v',  'QTY_V_WIND_COMPONENT',
      template_file = 'perfect_input.nc'
   /
   # only used if initial conditions file not specified in run
   &bgrid_cold_start_nml
      nlon = 60,
      nlat = 30,
      nlev = 5,
      equal_vert_spacing = .true.
   /
   # Values in hs_forcing_nml are described in Held and Suarez (1994)
   &hs_forcing_nml
      delh      =  60.,
      t_zero    = 315.,
      t_strat   = 200.,
      delv      =  10.,
      eps       =   0.,
      ka        = -40.,
      ks        =  -4.,
      kf        =  -1.,
      sigma_b   =  .7,
      do_conserve_energy = .false.
   /
   &bgrid_core_driver_nml
      damp_coeff_wind   = 0.10,
      damp_coeff_temp   = 0.10,
      damp_coeff_tracer = 0.10,
      advec_order_wind   = 4,
          advec_order_temp   = 2,
          advec_order_tracer = 2,
          num_sponge_levels = 1,
          sponge_coeff_wind   = 1.00,
          sponge_coeff_temp   = 1.00,
          sponge_coeff_tracer = 1.00,
          num_fill_pass = 2,
          decomp = 0,0,
          num_adjust_dt = 3,
          num_advec_dt  = 3,
          halo = 1,
          do_conserve_energy = .false.
   /
   &bgrid_integrals_nml
      file_name  = 'dynam_integral.out',
      time_units = 'days',
      output_interval = 1.00
   /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following values are specified in ``model_nml``.

+----------------------+--------------------+-------------------------------------------+
| Item                 | Type               | Description                               |
+======================+====================+===========================================+
| current_time(4)      | integer            | Specifies the initial time of the Bgrid   |
|                      |                    | model internal clock. The four integer    | 
|                      |                    | values are the day, hour, minutes, and    |
|                      |                    | seconds. The default version of the Bgrid |
|                      |                    | model has neither a diurnal or seasonal   |
|                      |                    | cycle, so these can all be set to 0, the  |
|                      |                    | default value.                            |
+----------------------+--------------------+-------------------------------------------+
| override             | logical            | If true, then the initial model date is   |
|                      |                    | taken from namelist entry current_time,   |
|                      |                    | even if an atmos_model.res file is found  |
|                      |                    | in directory INPUT. For most DART         |
|                      |                    | applications, atmospheric restart values  |
|                      |                    | are coming from DART files and no INPUT   |
|                      |                    | directory is used.                        |
+----------------------+--------------------+-------------------------------------------+
| dt_atmos             | integer            | Model timestep in seconds.                |
+----------------------+--------------------+-------------------------------------------+
| noise_sd             | real(r8)           | Standard deviation of random              |
|                      |                    | perturbations to the time tendency of     |
|                      |                    | temperature applied at each timestep.     |
|                      |                    | Each gridpoint value of the computed      |
|                      |                    | temperature tendency is multiplied by     |
|                      |                    | 1+N(0, noise_sd) before the updated       |
|                      |                    | values of temperature are computed.       |
+----------------------+--------------------+-------------------------------------------+
| dt_bias              | integer            | Allows a simple mechanism to simulate     |
|                      |                    | model error. If dt_bias is non-zero, the  |
|                      |                    | assimilation programs believe that each   |
|                      |                    | model advance changes the time by         |
|                      |                    | dt_bias. However, internally the bgrid    |
|                      |                    | model is moving things forward by         |
|                      |                    | dt_atmos. By running perfect_model_obs    |
|                      |                    | with one time step for the internal bgrid |
|                      |                    | clock (for instance dt_atmos = 3600,      |
|                      |                    | dt_bias = 3600), and filter with another  |
|                      |                    | (dt_atmos = 3000, and dt_bias = 3600)     |
|                      |                    | model error is simulated.                 |
+----------------------+--------------------+-------------------------------------------+
| state_variables(:,2) | character(len=129) | Strings that identify the bgrid_solo      |
|                      |                    | variables that should be part of the DART |
|                      |                    | state vector. The first column is the     | 
|                      |                    | netCDF variable name, the second column   |
|                      |                    | is the corresponding DART quantity.       |
+----------------------+--------------------+-------------------------------------------+
| template_file        | character(len=256) | This is the name of the file that         |
|                      |                    | specifies the resolution of the variables |
|                      |                    | DART uses to create the DART state        |
|                      |                    | vector. If *template_file = "null"* the   |
|                      |                    | *&bgrid_cold_start_nml* namelist          |
|                      |                    | variables are used to specify the         |
|                      |                    | resolution. The actual input filenames    |
|                      |                    | for *filter* and *perfect_model_obs* come |
|                      |                    | from their respective namelists.          |
|                      |                    | The resolutions in the file specified in  |
|                      |                    | *template_file* must match the            |
|                      |                    | resolutions of the variables in the input |
|                      |                    | filenames. To start an experiment with a  |
|                      |                    | new model resolution, set template_file   |
|                      |                    | to "null" and set the resolutions in      |
|                      |                    | bgrid_cold_start_nml.                     | 
+----------------------+--------------------+-------------------------------------------+

The following values are specified in ``bgrid_cold_start_nml``.

+------------------------+--------------------+-------------------------------------------+
| Item                   | Type               | Description                               |
+========================+====================+===========================================+
| nlon                   | integer            | The number of longitudes on the model     |
|                        |                    | grid.                                     |
+------------------------+--------------------+-------------------------------------------+
| nlat                   | integer            | The number of latitudes on the model      |
|                        |                    | grid.                                     |
+------------------------+--------------------+-------------------------------------------+
| nlev                   | integer            | The number of model levels.               |
+------------------------+--------------------+-------------------------------------------+
| equal_vertical_spacing | logical            | Model levels are equally spaced in        |
|                        |                    | pressure if true.                         |
+------------------------+--------------------+-------------------------------------------+

The Held-Suarez forcing details can be modified with the ``hs_forcing_nml``
namelist using the documentation in Held and Suarez (1994).

Model dynamics can be adjusted with the bgrid_core_driver_nml following the
documentation in the references and internal documentation in the bgrid code.

References
----------

.. [1] Anderson, J. L. and Coauthors, 2004: The new GFDL global atmosphere and
       land model AM2-LM2: Evaluation with prescribed SST simulations. *Journal
       of Climate*, **17**, 4641-4673. `doi:10.1175/JCLI-3223.1 <https://doi.org/10.1175/JCLI-3223.1>`_

.. [2] Held, I. M., and M. J. Suarez, 1994: A proposal for the intercomparison
       of the dynamical cores of atmospheric general circulation models,
       *Bulletin of the American Meteorological Society*, **75(10)**, 1825-1830.
       `doi:10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2 <https://doi.org/10.1175/1520-0477(1994)075\<1825:APFTIO\>2.0.CO;2>`_

.. [3] Anderson, J. L., Wyman, B., Zhang, S. & Hoar, T., 2005: Assimilation of
       surface pressure observations using an ensemble filter in an idealized
       global atmospheric prediction system, *Journal of the Atmospheric Sciences*,
       **62**, 2925-2938. `doi:10.1175/JAS3510.1 <https://doi.org/10.1175/JAS3510.1>`_

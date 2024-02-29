Radar observations
==================

Overview
--------

DART provides limited support for the conversion of radar observations to 
``obs_seq`` format. As an end goal, you want to assimilate radar observations
that:

* Have been quality controlled to remove non-meteorological scatterers and
  other artifacts
* Have horizontal resolution that has been reduced to approximately twice the
  expected horizontal grid spacing of your model. For example, if your model 
  has 3 km grid spacing, you should reduce your radar observations to every 6
  km, interpolated along the sweep plane.

Reflectivity observations are often partitioned into two types:

1. Regular reflectivity observations
2. Clear-air reflectivity observations where no radar echoes are observed.

Quality control is best done with raw data. You should have an ability to
perform quality control before converting your observations to obs_seq format.

Synthetic radar observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``create_obs_radar_sequence`` program generates one or more sets of 
synthetic `WSR-88D (NEXRAD) <http://en.wikipedia.org/wiki/WSR-88D>`__ radar
observations. It can generate reflectivity and/or doppler radial velocity
observations with clear-air or storm sweep patterns. These synthetic
observations can be used for testing your assimilation setup or for conducting
Observing System Simulation Experiments (OSSEs).

To build ``create_obs_radar_sequence``, change directory into the ``work`` 
subdirectory, ensure ``input.nml`` is configured properly and run the build
script:

.. code-block::

   cd work
   ./quickbuild.sh

Real radar observations
~~~~~~~~~~~~~~~~~~~~~~~

Once you have ensured that your data are quality controlled, use the 
`Observation Processing And Wind Synthesis (OPAWS) <http://code.google.com/p/opaws/>`__
utility convert your data to `obs_seq` format. The OPAWS utility reads specific
types of files as input, such as DORADE sweep files and NSF NCAR EOL Foray data.

OPAWS analyzes and grids data in either:

* two-dimensions (on the conical surface of each sweep), or
* three-dimensions (Cartesian).

If your raw data are not in such a format, additional utilities are available
for conversion such as the 
`RADX library which is part of the LIDAR/RADAR Open Software Environment <https://github.com/NCAR/lrose-core>`__.

Guidance for Weather Research and Forecasting (WRF) users
---------------------------------------------------------

If you intend to assimilate radar observations into WRF, you'll need to make
some code modifications to allow for forward operator calculations. For
reflectivity, most of the available microphysics schemes have built-in
capability to output reflectivity, assuming a 10 cm wavelength. If you are not
using an S-band radar, be aware that attenuation is not accounted for in the
built-in reflectivity operator.

For radial velocity, you will also need to generate a new diagnostic field:
terminal fall velocity. There is very limited support for fall velocity in WRF,
although it is partially supported in the Thompson microphysics scheme. 

.. note:: 

   You will still need to modify WRF code to get this diagnostic output to
   history files.

With these two fields available in your WRF history files, you can add them to
your DART `wrf_state_variables` list.

You should also use a special localization radius for radar observations, 
typically 12-24 km. If you leave range-folding in your radar observations, you
will need to build the special version of DART that unfolds the velocity
observations on-the-fly.

With all of those configurations in place, you will be ready to assimilate 
radar observations using WRF and DART.

For more information, see the WRF tests directory in
``DART/models/wrf/regression/Radar/`` for pointers to data to run a radar test
case.


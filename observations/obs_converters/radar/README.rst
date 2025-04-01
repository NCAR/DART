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
synthetic `WSR-88D (NEXRAD) <https://www.roc.noaa.gov/level-two-data-types.php>`__ 
radar observations. It can generate reflectivity and/or doppler radial velocity
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

The National Weather Service provides `NEXRAD (WSR-88D Type II) <https://www.roc.noaa.gov/level-two-data-types.php>`__
radar data for download. We recommend the use of the 
`Observation Processing And Wind Synthesis Python Tool (pyOPAWS) <https://github.com/louiswicker/pyOPAWS/>`__
to directly convert from NEXRAD to DART `obs_seq` format.  This tool was created by Dr.
Louis Wicker from NOAA, and also includes diagnostic scripting to visualize the data.  
Follow the `installation <https://github.com/louiswicker/pyOPAWS?tab=readme-ov-file#installation>`__ 
and `usage instructions <https://github.com/louiswicker/pyOPAWS?tab=readme-ov-file#usage-examples>`__ 
to compile and run the pyOPAWS software.  We also provide quickstart instructions to compile and run pyOPAWS
in the next section.  

If your radar data are not in NEXRAD format, additional conversion utilities are available
such as `LROSE (LIDAR/RADAR Open Software Environment) <https://github.com/NCAR/lrose-core>`__,
developed through a collaboration between the `Earth Observation Laboratory at NSF NCAR <http://lrose.net/>`__, 
and Colorado State University. The LROSE software, however, does not provide a direct converter
to DART `obs_seq` format.


A separate OPAWS utility is capable of converting DORADE sweep and NSF NCAR EOL Foray formatted data.
If you are working with these radar data formats please contact the DART
team directly (dart@ucar.edu).

pyOPAWS Quickstart Instructions for NSF NCAR's Derecho
------------------------------------------------------

In this example, pyOPAWS was run using a standard set of modules
as listed below. The pyOPAWS tool requires compilation and expects
GNU fortran (gcc) to be loaded.

::
 
  1) ncarenv/24.12 (S)      5) libfabric/1.15.2.0
  2) craype/2.7.31          6) cray-mpich/8.1.29
  3) gcc/12.4.0             7) hdf5/1.12.3
  4) ncarcompilers/1.0.0    8) netcdf/4.9.2

1.  Load the environment management system for python (conda) as: 

::

   module load conda

2. Navigate to your local pyOPAWS directory {pyOPAWS_DIR}, then activate the
   pyOPAWS python environment as:

::

   cd {pyOPAWS_DIR}
   conda env create -f pythnon4opaws.yml
   conda activate opaws_env

3. Confirm list of ``opaws_env`` packages:

::

  conda list --name opaws_env


  # packages in environment at conda-envs/opaws_env:
  #
  # Name                    Version                   Build  Channel
  _libgcc_mutex             0.1                 conda_forge    conda-forge
  _openmp_mutex             4.5                       2_gnu    conda-forge
  aiobotocore               2.21.1             pyhd8ed1ab_0    conda-forge
  aiohappyeyeballs          2.6.1              pyhd8ed1ab_0    conda-forge
  aiohttp                   3.11.13         py310h89163eb_0    conda-forge
  aioitertools              0.12.0             pyhd8ed1ab_1    conda-forge
  aiosignal                 1.3.2              pyhd8ed1ab_0    conda-forge
  alabaster                 1.0.0                    pypi_0    pypi
  ..
  .. and many more

4. To run the converter script (``opaws2D.py``) first compile the fortran cressman routine
   which assumes the GNU compiler was installed:

::

  python fcompile_cressman.py

Confirm ``cressman.cpython-310-x86_64-linux-gnu.so`` has been built in {pyOPAWS_DIR}.



5. (OPTIONAL) Follow the `example <https://github.com/louiswicker/pyOPAWS?tab=readme-ov-file#usage-examples>`__
   instructions to generate a diagnostic plot of radar data as:

::

   python lvl2_plot.py -f KOUN20110524_224700_V06 -q None -p 1 -u phase

View the file {pyOPAWS_DIR}/images/KOUN_20110524_224700_0.68DEG.png which contains
a multi-panel image of reflectivity, unfolded radial velocity, spectrum width,
Z-dr, RHO_H-V, PHI-DP. 
 

6. Follow the `example <https://github.com/louiswicker/pyOPAWS?tab=readme-ov-file#usage-examples>`__
   instructions to convert the NEXRAD file (KOUN20110524_224700_V06) to DART `obs_seq` format. 
   Note that the -w flag indicates the DART `obs_seq*.out` files will be generated.

::
  
  python opaws2d.py -f KOUN20110524_224700_V06 -q None -p 1 -w -u phase

Confirm the following files within the {pyOPAWS_DIR}/opaws_files/ were created:

::

  KOUN_20110524_224700_01.pdf
  obs_seq_KOUN_20110524_224700_VR.out
  obs_seq_KOUN_20110524_224700_RF.out
  obs_seq_KOUN_20110524_224700.nc

Note that two DART `obs_seq*.out` files have been generated that include
the radial velocity and reflectivity observations respectively. 

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


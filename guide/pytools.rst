.. index:: radiance, radiances, jedi

.. _pytools:

=========================
DART Python tools
=========================

Installation
------------

DART provides a set of Python tools that can be used to work with DART data and models.
These tools are located in the DART/pytools directory.

We recommend creating a `virtual environment <https://docs.python.org/3/library/venv.html>`__ 
to install DART Python tools. This helps to keep the DART Python tools and their dependencies 
isolated from other Python packages you may have installed.

.. code-block:: bash

    python3 -m venv py-dart
    source py-dart/bin/activate

You can then install the DART Python tools by navigating to the DART/pytools directory and running:

.. code-block:: bash

    cd DART/pytools
    pip install -r pyDART.txt


When you want to use the DART Python tools in a new session you will need to activate the virtual environment:

.. code-block:: bash

    source py-dart/bin/activate

To deactivate the virtual environment when you are done, simply run:

.. code-block:: text

    deactivate

Available Tool Packages
-----------------------

.. _ioda2obsq:

pyjedi
^^^^^^

The `pyjedi` package provides utilities for DART JEDI interoperability.
This is a new package that is under development and is not fully functional yet.
At this time, a converter ``ioda2obsq`` is available that can transform a JEDI IODA observation file into a DART observation sequence file.
This converter can be run as follows:

.. code-block:: bash

    ioda2obsq <yamlConfigFile> <inputIodaFile> <outputObsqFile>

where ``<yamlConfigFile>`` is a YAML configuration file that specifies the mapping between JEDI IODA observation
variables and DART observation variables, ``<inputIodaFile>`` is the input JEDI IODA observation file,
and ``<outputObsqFile>`` is the output DART observation sequence file.

To view the ioda2obsq help message, run:

.. code-block:: bash

  ioda2obsq -h


The YAML configuration file contains the mapping between JEDI IODA observation variables and DART observation types along
with the specification of the desired vertical coordinate variable.
Here is an example configuration file for radiosondes:

.. code-block:: yaml

    ---

    ioda to obsq converter:
      observation variables:
        - name: airTemperature
          type: RADIOSONDE_TEMPERATURE
        - name: specificHumidity
          type: RADIOSONDE_SPECIFIC_HUMIDITY
        - name: windEastward
          type: RADIOSONDE_U_WIND_COMPONENT
        - name: windNorthward
          type: RADIOSONDE_V_WIND_COMPONENT

      observation category:
        name: conventional
        vertical coordinate:
          name: MetaData/pressure
          units: "pressure (Pa)"

The ``observation category:`` section is how the tool discerns between different categories of observation types.
Currently, the tool recognizes "conventional" (radiosonde, aircraft, etc.) and "radiance" (amsu-a, goes, etc.) with the ``name:`` specification.
The idea with the observation category is to group instruments where the data are organized in a similar fashion.
For example, in the future, perhaps "radar" and "gpsro" could be added to the list of categories.

The ``name:`` entries for variables are the names of the variables in the JEDI IODA observation file, and the ``type:`` entries
are the corresponding DART observation types.
The ``vertical coordinate:`` section specifies the variable that will be used as the vertical coordinate in
the DART observation sequence file.
Note that in the ``vertical coordinate:`` section, the ``name:`` entry must specify the full netcdf path to the variable
in the JEDI IODA observation file, and the ``units:`` entry specifies the units to be used in the DART observation sequence file.

Here is an example radiance obs type (AMSU-A):

.. code-block:: yaml

  ---

  ioda to obsq converter:
    observation variables:
      - name: brightnessTemperature
        type: NOAA_19_AMSUA_TB

    observation category:
      name: radiance
      channel numbers: 1, 2, 3, 12-15
      vertical coordinate:
        units: "pressure (Pa)"
        data value: 35000.0
      metadata:
        sensor key: NOAA_19_AMSUA
        rttov sensor db: DART/observations/forward_operators/rttov_sensor_db.csv
        sat az variable: MetaData/sensorAzimuthAngle
        sat ze variable: MetaData/sensorZenithAngle

Some of the entries are similar to those in the radiosonde (conventional) category.
Note that in the radiance category, the ``vertical coordinate:`` spec needs a ``units:`` (as before) and a ``data value:``. The 
``data value:`` is used to fill in the vertical coordinate variable for all of the observations since radiance observations 
do not have a vertical coordinate associated with them. Currently, the ``data value:`` is repeated on all of the observations.

Two new entries under the ``observation category:`` spec have been introduced.
The ``channel numbers:`` spec allows the user to select a subset of channels to convert.
The value for ``channel numbers:`` is a comma separated list where each entry can be an integer or a range of integers.

The ``metadata:`` spec allows the user to configure what gets placed in the obs sequence file's metadata section.
The ``sensor key:`` and ``rttov sensor db:`` spec go hand-in-hand, and describe how to pull out id numbers for the platform, satellite and sensor.
``rttov_sensor_db.csv`` is a DART file that contains the mapping between sensor keys and their corresponding
platform, satellite and sensor id numbers. Refer to :ref:`obs_def_rttov_mod` for more detail.
The ``sat az variable:`` and ``sat ze variable:`` are the JEDI IODA names for the satellite azimuth and zenith angles respectively.
Another piece of information that comes from the ``rttov_sensor_db:`` file is the spectral band (eg. infrared, microwave).
Currently, only infrared (``ir``) and microwave (``mw``) are recognized.
The purpose of the ``metadata:`` spec is to have flexibility with the configuration for infrared, microwave, visible, etc. instruments.

.. note::

    The ``ioda2obsq`` tool is under active development and has limited functionality at this time.
    It is expected that more features will be added soon, including support for additional satellite radiance observation types.
    In its current state, it is primarily intended for use with radiosonde and similar conventional observation types, as well as infrared and microwave radiance instruments.

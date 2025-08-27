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

pyjedi
^^^^^^

The `pyjedi` package provides utilities for DART JEDI interoperability.
This is a new package that is under development and is not fully functional yet.
At this time, a new converter is available that can transform a JEDI IODA observation file into a DART observation sequence file.
This converter can be run as follows:

.. code-block:: bash

    ioda2obsq <yamlConfigFile> <inputIodaFile> <outputObsqFile>

where ``<yamlConfigFile>`` is a YAML configuration file that specifies the mapping between JEDI IODA observation
variables and DART observation variables, ``<inputIodaFile>`` is the input JEDI IODA observation file,
and ``<outputObsqFile>`` is the output DART observation sequence file.

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

      vertical coordinate:
        name: MetaData/pressure
        units: "pressure (Pa)"

The `name:` entries are the names of the variables in the JEDI IODA observation file, and the `type:` entries
are the corresponding DART observation types.
The `vertical coordinate` section specifies the variable that will be used as the vertical coordinate in
the DART observation sequence file.
Note that in the `vertical coordinate` section, the `name:` entry must specify the full netcdf path to the variabile
in the JEDI IODA observation file, and the `units:` entry specifies the units to be used in the DART observation sequence file.

.. note::

    The `ioda2obsq` tool is under active development and has limited functionality at this time.
    It is expected that more features will be added soon, including support for satellite radiance observation types.
    In its current state, it is primarily intended for use with radiosonde and similar conventional observation types.
PROGRAM ``cice_to_obs``
=======================

Overview
--------

Sea ice percentage observations to DART converter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This converter reads the binary sea ice observations from the snow and ice data center files and outputs DART obs_seq
format files. It will loop over multiple days inside a single run of the converter program.

Data sources
------------

The `National Snow and Ice Data Center <http://nsidc.org/>`__ supplies the data files read by this converter. (I believe
it is `this format? <http://nsidc.org/data/NSIDC-0051>`__)

Programs
--------

The ``cice_to_obs.f90`` file is the source for the main converter program. More documentation is in the source code file
especially around where the namelist variables are declared.

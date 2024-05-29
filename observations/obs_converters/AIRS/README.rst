AIRS and AMSU
=============

The AIRS directory contains both the AIRS and AMSU-A observation converters.
AIRS is the Atmospheric Infrared Sounder (AIRS) Level 2 observations.
AMSU-A is the Advanced Microwave Sounding Unit (AMSU-A) L1B Brightness Temperatures.

- :doc:`./convert_airs_L2` is used to convert AIRS temperature and 
  specific humidity vertical profile observations.
- :doc:`./convert_amsu_L1` is used to convert AMSU-A radiances (brightness temperature) 
  observations.

Both converters are in the AIRS directory because of the complicated history
of the data used to create the AIRS L2 product (which includes some AMSU observations).
Since both datasets are HDF - it was believed that some of the routines could be
used by both converters. Alas, that has not proven to be the case.


Dependencies
------------

Both ``convert_airs_L2`` and ``convert_amsu_L1`` require the HDF-EOS2 libraries,
which, in turn, requires HDF4. HDF4 is available on Derecho using the ``module load hdf``
command.

The ``convert_amsu_L1`` script also requires the RTTOV libraries.

The following mkmf.templates for gfortran and intel compilers respectively, 
are available in DART/build_templates, and they have been used to compile 
the AIRS and AMSUA observation converters on Derecho. They include the 
proper library paths for both HDF-EOS2 and RTTOV. The HDF-EOS2 library 
required a patch to work with DART observation converters.
For details on the patch see `issue #590 <https://github.com/NCAR/DART/issues/590>`_.

.. code :: text
 
    mkmf.template.AIRS.gfortran
    mkmf.template.AIRS.intel

In addition to gfortran and intel compiled hdf-eos as mentioned above, 
Derecho also includes nvhpc and cray compiled hdf-eos libraries, with 
the paths provided below.

.. code:: text

    /glade/campaign/cisl/dares/libraries/hdf-eos_intel
    /glade/campaign/cisl/dares/libraries/hdf-eos_gfortran
    /glade/campaign/cisl/dares/libraries/hdf-eos_nvhpc
    /glade/campaign/cisl/dares/libraries/hdf-eos_cray


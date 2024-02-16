AIRS and AMSU
=============

The AIRS directory contains both the AIRS and AMSU-A observation converters.
AIRS is the Atmospheric Infrared Sounder (AIRS) Level 2 observations.
AMSU-A is the Advanced Microwave Sounding Unit (AMSU-A) L1B Brightness Temperatures.

- :doc:`./convert_airs_L2` for converting AIRS temperature and moisture retrievals.
- :doc:`./convert_amsu_L1` for converting AMSU-A radiances.

Both converters are in the AIRS directory because of the complicated history
of the data used to create the AIRS L2 product (which includes some AMSU observations).
Since both datasets are HDF - it was believed that some of the routines could be
used by both converters. Alas, that has not proven to be the case.


Dependencies
------------

Both ``convert_airs_L2`` and ``convert_amsu_L1`` require the HDF-EOS2 libraries.

The HDF-EOS2 library required a patch to work with DART observation converters.
For details on the patch see `issue #590 <https://github.com/NCAR/DART/issues/590>`_.
The patched version of the library is available on Derecho. The _intel, _gfortran,
_nvhpc, _cray indicates which compiler was used to build the library. Select an
appropriate mkmf.template for the compiler you are using.

.. code:: text

    /glade/campaign/cisl/dares/libraries/hdf-eos_intel
    /glade/campaign/cisl/dares/libraries/hdf-eos_gfortran
    /glade/campaign/cisl/dares/libraries/hdf-eos_nvhpc
    /glade/campaign/cisl/dares/libraries/hdf-eos_cray

hdf-eos requires HDF4. HDF4 is available on Derecho with ``module load hdf4``.

``convert_amsu_L1`` requires the RTTOV libraries.

The following mkmf.templates, which are available in DART/build_templates have been used to compile the AIRS and AMSUA
observation converters on Derecho.

.. code :: text
 
    mkmf.template.AIRS.gfortran
    mkmf.template.AIRS.intel




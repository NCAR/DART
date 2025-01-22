QuikSCAT SeaWinds Data
======================

Overview
--------

NASA's QuikSCAT mission is described in
`Quick Scatteromoeter <https://podaac.jpl.nasa.gov/QuikSCAT>`_. "QuikSCAT"
refers to the satellite, "SeaWinds" refers to the instrument that provides near-surface wind speeds and directions over
large bodies of water. QuikSCAT has an orbit of about 100 minutes, and the SeaWinds microwave radar covers a swath under
the satellite. The swath is comprised of successive scans (or rows) and each scan has many wind-vector-cells (WVCs). For
the purpose of this document, we will focus only the **Level 2B** product at 25km resolution. If you go to the official
JPL data distribution site `podaac.jpl.nasa.gov <http://podaac.jpl.nasa.gov/DATA_CATALOG/quikscatinfo.html>`_, we are using the product labelled
**L2B OWV 25km Swath**. Each orbit consists of (potentially) 76 WVCs in each of 1624 rows or scans. The azimuthal
diversity of the radar returns affects the error characteristics of the retrieved wind speeds and directions, as does
rain, interference of land in the radar footprint, and very low wind speeds. Hence, not all wind retrievals are created
equal.

The algorithm that converts the 'sigma naughts' (the measure of radar backscatter) into wind speeds and directions has
multiple solutions. Each candidate solution is called an 'ambiguity', and there are several ways of choosing 'the best'
ambiguity. At present, the routine to convert the original L2B data files (one per
orbit) in HDF format into the DART observation sequence file makes several assumptions:

#. All retrievals are labelled with a 10m height, in accordance with the retrieval algorithm.
#. Only the highest-ranked (by the MLE method) solution is desired.
#. Only the WVCs with a wvc_quality_flag of **zero** are desired.
#. The mission specification of a wind speed rms error of 2 ms (for winds less than 20 m/s) and 10% for windspeeds
   between 20 and 30 m/s can be extended to all winds with a qc flag of zero.
#. The mission specification of an error in direction of 20 degrees rms is applicable to all retrieved directions.
#. All retrievals with wind speeds less than 1.0 are not used.
#. The above error characterstics can be simplified when deriving the horizontal wind components (i.e. U,V). **Note :**
   this may or may not be a good assumption, and efforts to assimilate the speed and direction directly are under way.

Data sources
------------

The NASA Jet Propulsion Laboratory (JPL) `data repository <https://podaac.jpl.nasa.gov/>`_ has a
collection of animations and data sets from this instrument. In keeping with NASA tradition, these data are in HDF
format (specifically, HDF4), so if you want to read these files directly, you will need to install the 
`HDF4 libraries <https://portal.hdfgroup.org/display/support/Download+HDF4>`_.

If you go to the official JPL data distribution site http://podaac.jpl.nasa.gov/DATA_CATALOG/quikscatinfo.html, we are
using the product labelled **L2B OWV 25km Swath**. They are organized in folders by day, with each orbit (each
revolution) in one compressed file. There are 14 revolutions per day. The conversion to DART observation sequence format
is done on each revolution, multiple revolutions may be combined 'after the fact' by any ``obs_sequence_tool`` in the
``work`` directory of any model.

Programs
--------

There are several programs that are distributed from the JPL www-site,
ftp://podaac.jpl.nasa.gov/pub/ocean_wind/quikscat/L2B/sw/; we modified the Fortran file
`read_qscat2b.f <ftp://podaac.jpl.nasa.gov/pub/ocean_wind/quikscat/L2B/sw/FORTRAN/read_qscat2b.f>`__ 
to be a subroutine for use with DART. For reference, the original ``read_qscat2b.f`` and ``Makefile``
are included in ``DART/observations/quikscat`` directory.


convert_L2b.f90
~~~~~~~~~~~~~~~

``convert_L2b`` converts the HDF files distributed by JPL to an obs_sequence file.
To build ``convert_l2b`` using ``quickbuild.sh`` you will need the HDF4 library.
HDF4 is available on the NSF NCAR machine Derecho: ``module load hdf``.

.. warning::

  To avoid conflicts with netCDF library required by DART, we recommend building HDF4 *without* 
  the HDF4 versions of the NetCDF API. 

After successfully building HDF, add the appropriate library flags to your mkmf.template file,
or for Derecho users, use the already available files at 
*DART/build_templates/mkmf.template.quikscat.intel* or 
*DART/build_templates/mkmf.template.quikscat.gfortran*. Below is a snippet from an
mkmf.template file used to link to both NetCDF and HDF4.   

.. code:: text

   NETCDF = /glade/u/apps/derecho/23.06/spack/opt/spack/netcdf/4.9.2/oneapi/2023.0.0/iijr
   HDF = /glade/u/apps/derecho/23.06/spack/opt/spack/hdf/4.2.15/oneapi/2023.0.0/yo2r
   
   INCS = -I$(NETCDF)/include -I$(HDF)/include
   LIBS = -L$(NETCDF)/lib -lnetcdff -lnetcdf -L$(HDF)/lib -lmfhdf -ljpeg -lz -ltirpc -lsz -ldf
   FFLAGS  = -O -assume buffered_io $(INCS)
   LDFLAGS = $(FFLAGS) $(LIBS)


There are a lot of observations in every QuikSCAT orbit. Consequently, the observation sequence files are pretty large -
particularly if you use the ASCII format. Using the binary format (i.e. *obs_sequence_nml:write_binary_obs_sequence =
.true.*) will result in observation sequence files that are about *half* the size of the ASCII format.

Since there are about 14 QuikSCAT orbits per day, it may be useful to convert individual orbits to an observation
sequence file and then concatenate multiple observation sequence files into one file per day. This can be
accomplished with the :ref:`obs_sequence_tool<obs sequence tool>` program. To build the ``obs_sequence_tool``, 
add ``obs_sequence_tool`` to the list of programs in ``quickbuid.sh``.


Namelist
--------

This namelist is read from the file ``input.nml``. We adhere to the F90 standard of starting a namelist with an
ampersand '&' and terminating with a slash '/' for all our namelist input. Character strings that contain a '/' must be
enclosed in quotes to prevent them from prematurely terminating the namelist. The following values are the defaults for
these namelist items.

::

   &convert_L2b_nml
      l2b_file = '',
      datadir = '.',
      outputdir = '.',
      lon1 = 0.0, 
      lon2 = 360.0, 
      lat1 = -90.0, 
      lat2 = 90.0,
      along_track_thin = 0,
      cross_track_thin = 0
    /

| 

.. container::

   It is possible to restrict the output observation sequence to contain data from a region of interest throught the use
   of the namelist parameters. If you need a region that spans the Prime Meridian lon1 can be a larger number than lon2,
   for example, a region from 300 E to 40 E and 60 S to 30 S (some of the South Atlantic), would be *lon1 = 300, lon2 =
   40, lat1 = -60, lat2 = -30*.

   +------------------+--------------------+----------------------------------------------------------------------------+
   | Contents         | Type               | Description                                                                |
   +==================+====================+============================================================================+
   | l2b_file         | character(len=128) | name of the HDF file to read - NOT including the directory, e.g.           |
   |                  |                    | QS_S2B44444.20080021548                                                    |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | datadir          | character(len=128) | the directory containing the HDF files                                     |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | outputdir        | character(len=128) | the directory for the output observation sequence files.                   |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | lon1             | real(r4)           | the West-most longitude of interest in degrees. [0.0, 360]                 |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | lon2             | real(r4)           | the East-most longitude of interest in degrees. [0.0, 360]                 |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | lat1             | real(r4)           | the South-most latitude of interest in degrees. [-90.0, 90.0]              |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | lat2             | real(r8)           | the North-most latitude of interest in degrees. [-90.0, 90.0]              |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | along_track_thin | integer            | provides ability to thin the data by keeping only every Nth row. e.g. 3 == |
   |                  |                    | keep every 3rd row.                                                        |
   +------------------+--------------------+----------------------------------------------------------------------------+
   | cross_track_thin | integer            | provides ability to thin the data by keeping only every Nth wind vector    |
   |                  |                    | cell in a particular row. e.g. 5 == keep every 5th cell.                   |
   +------------------+--------------------+----------------------------------------------------------------------------+

|

Future Plans
~~~~~~~~~~~~

1. There is one bit of error-checking that did not survive the conversion from F77 to F90. I need to restore the check that the HDF file being read is a 'Level 2B' product.
2. There is a lot of error-checking that is not being done. I need to bulletproof the code more.
3. We need namelist options to select something other than the highest-ranked ambiguity.
4. We need namelist options to select more QC flags - not just the ones with the 'perfect' QC value of 0
5. Add an option to leave the observations as speed and direction instead of converting them to U,V components. This is a natural implementation of the instrument error characteristics. However, it would require writing a specialized forward operator in order to assimilate obs of this type (speed, direction), and there is still a numerical problem with trying to do the statistics required during the assimilation of a cyclic direction value.
 

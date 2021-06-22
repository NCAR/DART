#######################
Locating netCDF library
#######################

DART uses the `netCDF <https://www.unidata.ucar.edu/software/netcdf/>`__
self-describing data format for storing the results of assimilation
experiments. These files have the extension *.nc* and can be read by a
number of standard data analysis tools. In particular, DART also makes
use of the F90 netCDF interface which is available through the
``netcdf.mod`` and ``typesizes.mod`` modules and the ``libnetcdf``
library. Depending on the version, the ``libnetcdff`` library is also
often required.

If the netCDF library does not exist on your system, you must build it
(as well as the F90 interface modules).

.. warning::

   You must build netCDF with the same compiler (including version) you plan to
   use for compiling DART. In practice this means that even if you have a netCDF
   distribution on your system, you may need to recompile netCDF in a separate
   location to match the compiler you will use for DART. The library and
   instructions for building the library or installing from a package manager
   may be found at the
   `netCDF home page <https://www.unidata.ucar.edu/software/netcdf/>`_.

.. important::

   The normal location for the netCDF Fortran modules and libraries would be in
   the ``include`` and ``lib`` subdirectories of the netCDF installation.
   However, different compilers or package managers sometimes place the modules
   and/or libraries into non-standard locations. It is required that both
   modules and the libraries be present.

.. note::

   The location of the netCDF library, ``libnetcdf.a``, and the locations of
   both ``netcdf.mod`` and ``typesizes.mod`` will be needed later. Depending on
   the version of netCDF and the build options selected, the Fortran interface
   routines may be in a separate library named ``libnetcdff.a`` (note the two
   F’s). In this case both libraries are required to build executables.
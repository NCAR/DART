##################
Installing NetCDF 
##################

DART requires the `NetCDF libraries <https://www.unidata.ucar.edu/software/netcdf/>`__.
NetCDF is a set of software libraries and machine-independent data formats 
that support the creation, access, and sharing of array-oriented scientific data.

DART Requirements:

- NetCDF-C library (version ≥ 4.0)
- NetCDF-Fortran library (compatible with your compiler and NetCDF-C version)

Instructions for install netcdf can be found at the
`NetCDF home page <https://www.unidata.ucar.edu/software/netcdf/>`_.
Note, most Supercomputers and HPC clusters have a NetCDF module available,
and NetCDF is available from many package managers (e.g., apt, yum, macports, conda).

Ensure that both NetCDF libraries (C and Fortran) are compiled with support for 
your preferred compiler (e.g., gfortran, ifx).

Once you have installed netCDF, you will need to tell DART where to find the
NetCDF headers and libraries.

.. code-block:: bash

   nf-config --all  

will show you the NetCDF install location (prefix) and the paths to the NetCDF 
Fortran headers and libraries.
# template for the Intel Fortran Compiler Version 7.1 on machine "dart"
# typical use with mkmf
# mkmf -t template.ifc -c"-Duse_libMPI -Duse_netCDF" path_names /usr/local/include
#
# LIBS needs to be customized for your site
FC = ifc
LD = ifc
FFLAGS = -i4 -r8 -fpp -O2 -I/opt/Intel/include
LIBS = -L/usr/local/lib -lPEPCF90 -L/opt/Intel/lib -ludunits -lnetcdf -L/usr/mpi/lib
LDFLAGS = $(LIBS)

#
# <next four lines automatically updated by CVS, do not edit>
# $Source$
# $Revision$
# $Date$
# $Author$
# $Id$
#

# These commented-out scenarios are not guaranteed in any way ...
#
#LIBS		= -L/usr/local/lib32/r4i4 -lnetcdf -lm
#INCL		= -I/usr/local/include -I./
#FC		= xlf
#FCFLAGS	= -g -C -qfree=f90 -qarch=auto -qzerosize -qnosave -qmaxmem=-1 -qspill=20000

#LIBS		= -L/usr/local/netcdf/lib -lnetcdf -lm
#INCL		= -I/usr/local/netcdf/include -I./
#FC		= f90
#FCFLAGS	= -C -free -fpe1 -g

# LIBS		= -L/usr/local/netcdf-pgi/lib -lnetcdf -lm
# INCL		= -I/usr/local/netcdf-pgi/include -I./
# FC		= pgf90
# FCFLAGS	= -C -Mfree -Ktrap=fp -tp p6  -byteswapio

LIBS		= -L/usr/local/netcdf-3.5.1-beta11-pgi-hpf-cc-5.0-2a/lib/ -lnetcdf
INCL		= -I/usr/local/netcdf-3.5.1-beta11-pgi-hpf-cc-5.0-2a/include/
FC		= pgf90
FCFLAGS		= -pc 64 -v -C -Ktrap=fp

OBJS		= update_wrf_bc.o		\
		module_netcdf_interface.o	\
		module_couple_uv.o

update_wrf_bc:  $(OBJS)
		$(FC) -o $@ $(FCFLAGS) $(OBJS) $(LIBS)

update_wrf_bc.o: module_netcdf_interface.o       \
               module_couple_uv.o              \
               update_wrf_bc.f90
		$(FC) -c $(FCFLAGS) update_wrf_bc.f90 $(INCL)

module_netcdf_interface.o: module_netcdf_interface.f90
		$(FC) -c $(FCFLAGS) module_netcdf_interface.f90 $(INCL)

module_couple_uv.o: module_couple_uv.f90
		$(FC) -c $(FCFLAGS) module_couple_uv.f90 $(INCL)

clean:	
	rm -f *.o *.mod update_wrf_bc

# compaq flags

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

 LIBS		= -L/usr/local/netcdf-intel/lib -lnetcdf -lm
 INCL		= -I/usr/local/netcdf-intel/include -I./
 FC		= ifc
 FCFLAGS	= -FR -cm -w -I.

OBJS		= \
		update_wrf_bc.o			\
		module_netcdf_interface.o	\
		module_couple_uv.o

update_wrf_bc:	$(OBJS)
		$(FC) -o $@ $(FCFLAGS) $(OBJS) $(LIBS)

update_wrf_bc.o: \
		module_netcdf_interface.o	\
		module_couple_uv.o		\
		update_wrf_bc.F
		$(FC) -c $(FCFLAGS) update_wrf_bc.F $(INCL)

module_netcdf_interface.o: \
		module_netcdf_interface.F
		$(FC) -c $(FCFLAGS) module_netcdf_interface.F $(INCL)

module_couple_uv.o: \
		module_couple_uv.F
		$(FC) -c $(FCFLAGS) module_couple_uv.F $(INCL)

clean:	
	rm -f *.o *.mod update_wrf_bc


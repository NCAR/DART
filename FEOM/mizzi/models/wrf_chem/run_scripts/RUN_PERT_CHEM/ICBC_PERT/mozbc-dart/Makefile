.SUFFIXES: .o .f90

F90    = $(FC)
LIBS   = -L$(NETCDF_DIR)/lib $(NETCDF_AR_FILES)
INCLUDE_MODULES = -I$(NETCDF_DIR)/include

EXEC = mozbc 

OBJS = mo_calendar.o\
       mo_utils.o\
       mo_mozart_lib.o \
       mo_wrfchem_lib.o \
       main_bc_wrfchem.o

##dependencies
#$(OBJECTS) : makefile

${EXEC} :       ${OBJS}
		${F90} -o $@ ${FFLAGS} ${OBJS} ${LIBS} 

.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

cleanup:
		rm -f ${OBJS} *.mod

clean:
		rm -f core ${EXEC} ${OBJS} *.mod

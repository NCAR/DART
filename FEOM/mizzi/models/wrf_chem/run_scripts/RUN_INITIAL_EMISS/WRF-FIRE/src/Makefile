.SUFFIXES: .o .f90

F90    = $(FC)
LIBS   = -L$(NETCDF_DIR)/lib $(NETCDF_AR_FILES)
INCLUDE_MODULES = -I$(NETCDF_DIR)/include

EXEC = fire_emis 

OBJS = attr_types.o \
       fire_utils.o\
       wrf_utils.o\
       glb_utils.o\
       fire_srf_types.o \
       fire_file.o \
       fire_emis.o

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

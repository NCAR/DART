.SUFFIXES: .o .f90

F90    = $(FC)
LIBS   = -L$(NETCDF_DIR)/lib $(AR_FILES)
INCLUDE_MODULES = -I$(NETCDF_DIR)/include

ifeq ($(FC),lf95)
  FFLAGS = --g
else
  FFLAGS = -g
endif

ifeq ($(FC),pgf90)
  FFLAGS += -Mnosave -Ktrap=fp
endif

ifeq ($(FC),gfortran)
  FFLAGS += -ffree-line-length-none
endif

ifeq ($(FC),lf95)
  FFLAGS += --nap --nchk --pca --nsav --trace --trap --wide
endif

ifeq ($(FC),ifort)
  FFLAGS += -fpe0 -ftrapuv
endif

EXEC = mozbc 

OBJS = mo_calendar.o\
       mo_utils.o\
       mo_mozart_lib.o \
       mo_wrfchem_lib.o \
       main_bc_wrfchem.o

##dependencies
#$(OBJECTS) : makefile

${EXEC} :       ${OBJS}
		${F90} -o $@ ${OBJS} ${LIBS} 

.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

cleanup:
		rm -f ${OBJS} *.mod

clean:
		rm -f core ${EXEC} ${OBJS} *.mod

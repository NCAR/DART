.SUFFIXES: .o .f90

F90    = $(FC)
LIBS   = -L$(NETCDF_DIR)/lib $(AR_LIBS)
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

OBJS_BIO_EMISS = misc_definitions_module.o\
       constants_module.o\
       bio_types.o\
       area_mapper.o\
       bio_emiss.o

OBJS_SURF_XFORM = surfdata_xform.o

OBJS_MEGAN_XFORM = megan_xform.o

.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

##dependencies

megan_bio_emiss: ${OBJS_BIO_EMISS}
		 ${F90} -o $@ ${OBJS_BIO_EMISS} ${LIBS} 

surfdata_xform: ${OBJS_SURF_XFORM}
		${F90} -o $@ ${OBJS_SURF_XFORM} ${LIBS} 

megan_xform:    ${OBJS_MEGAN_XFORM}
		${F90} -o $@ ${OBJS_MEGAN_XFORM} ${LIBS} 

cleanup:
		rm -f *.o *.mod

clean:
		rm -f core *.o *.mod

.SUFFIXES: .o .f90

F90    = $(FC)
LIBS   = -L$(NETCDF_DIR)/lib -lnetcdf
INCLUDE_MODULES = -I$(NETCDF_DIR)/include
FFLAGS = -g

OBJS_WESELY = mo_calendar.o\
              make_wesely.o

OBJS_EXO_COLDENS = mo_calendar.o\
                   make_exo_coldens.o

.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

##dependencies

exo_coldens: ${OBJS_EXO_COLDENS}
		 ${F90} -o $@ ${FFLAGS} ${OBJS_EXO_COLDENS} ${LIBS} 

wesely: ${OBJS_WESELY}
		${F90} -o $@ ${FFLAGS} ${OBJS_WESELY} ${LIBS} 

cleanup:
		rm -f *.o *.mod

clean:
		rm -f core *.o *.mod

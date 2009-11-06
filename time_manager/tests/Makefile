

# get the flag settings from the mkmf template directly

include ../../mkmf/mkmf.template

OBJS = types_mod.o utilities_mod.o time_manager_mod.o null_mpi_utilities_mod.o
MODS = types_mod.mod utilities_mod.mod time_manager_mod.mod mpi_utilities_mod.mod
EXES = timetest
LOGS = dart_log.out dart_log.nml

all: $(EXES)

check: $(EXES)
	./timetest

clean:
	rm -fr $(OBJS) $(MODS) $(EXES) timetest.o $(LOGS)

timetest: $(OBJS) timetest.o
	ifort -o timetest timetest.o $(OBJS) $(LDFLAGS)

types_mod.o: ../../common/types_mod.f90
	ifort -c ../../common/types_mod.f90 $(FFLAGS)

null_mpi_utilities_mod.o: ../../mpi_utilities/null_mpi_utilities_mod.f90 types_mod.o
	ifort -c ../../mpi_utilities/null_mpi_utilities_mod.f90 $(FFLAGS)

utilities_mod.o: ../../utilities/utilities_mod.f90 types_mod.o
	ifort -c ../../utilities/utilities_mod.f90 $(FFLAGS)

time_manager_mod.o: ../time_manager_mod.f90 utilities_mod.o types_mod.o
	ifort -c ../time_manager_mod.f90 $(FFLAGS)

timetest.o: timetest.f90 time_manager_mod.o
	ifort -c timetest.f90 $(FFLAGS)


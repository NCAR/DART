# Makefile created by mkmf $Id: mkmf 4255 2010-02-05 00:39:42Z thoar $ 

include ../../mkmf/mkmf.template

CPPDEFS = -Duse_netCDF

.DEFAULT:
	-touch $@
all: nml_test
time_manager_mod.o: ../../time_manager/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../time_manager/time_manager_mod.f90
null_mpi_utilities_mod.o: ../../mpi_utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../mpi_utilities/null_mpi_utilities_mod.f90
nml_test.o: ../../utilities/nml_test.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../utilities/nml_test.f90
utilities_mod.o: ../../utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../utilities/utilities_mod.f90
types_mod.o: ../../common/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../common/types_mod.f90
./time_manager_mod.f90: ../../time_manager/time_manager_mod.f90
	cp ../../time_manager/time_manager_mod.f90 .
./null_mpi_utilities_mod.f90: ../../mpi_utilities/null_mpi_utilities_mod.f90
	cp ../../mpi_utilities/null_mpi_utilities_mod.f90 .
./utilities_mod.f90: ../../utilities/utilities_mod.f90
	cp ../../utilities/utilities_mod.f90 .
./types_mod.f90: ../../common/types_mod.f90
	cp ../../common/types_mod.f90 .
./nml_test.f90: ../../utilities/nml_test.f90
	cp ../../utilities/nml_test.f90 .
SRC = ../../time_manager/time_manager_mod.f90 ../../mpi_utilities/null_mpi_utilities_mod.f90 ../../utilities/nml_test.f90 ../../utilities/utilities_mod.f90 ../../common/types_mod.f90
OBJ = time_manager_mod.o null_mpi_utilities_mod.o nml_test.o utilities_mod.o types_mod.o
OFF = ../../time_manager/time_manager_mod.f90 ../../mpi_utilities/null_mpi_utilities_mod.f90 ../../utilities/utilities_mod.f90 ../../common/types_mod.f90 ../../utilities/nml_test.f90
clean: neat
	-rm -f .cppdefs $(OBJ) nml_test
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
nml_test: $(OBJ)
	$(LD) $(OBJ) -o nml_test  $(LDFLAGS)

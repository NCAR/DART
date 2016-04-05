# Makefile created by mkmf $Id$ 

include ../../../mkmf/mkmf.template


.DEFAULT:
	-touch $@
all: fixsys preprocess
fixsys: 
	@ (cd ../../../mpi_utilities; ./fixsystem $(FC) )
time_manager_mod.o: ../../../time_manager/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../time_manager/time_manager_mod.f90
null_mpi_utilities_mod.o: ../../../mpi_utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../mpi_utilities/null_mpi_utilities_mod.f90
preprocess.o: ../../../preprocess/preprocess.f90 utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../preprocess/preprocess.f90
utilities_mod.o: ../../../utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../utilities/utilities_mod.f90
types_mod.o: ../../../common/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../common/types_mod.f90
./utilities_mod.f90: ../../../utilities/utilities_mod.f90
	cp ../../../utilities/utilities_mod.f90 .
./types_mod.f90: ../../../common/types_mod.f90
	cp ../../../common/types_mod.f90 .
./preprocess.f90: ../../../preprocess/preprocess.f90
	cp ../../../preprocess/preprocess.f90 .
./null_mpi_utilities_mod.f90: ../../../mpi_utilities/null_mpi_utilities_mod.f90
	cp ../../../mpi_utilities/null_mpi_utilities_mod.f90 .
./time_manager_mod.f90: ../../../time_manager/time_manager_mod.f90
	cp ../../../time_manager/time_manager_mod.f90 .
SRC = ../../../time_manager/time_manager_mod.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../preprocess/preprocess.f90 ../../../utilities/utilities_mod.f90 ../../../common/types_mod.f90
OBJ = time_manager_mod.o null_mpi_utilities_mod.o preprocess.o utilities_mod.o types_mod.o
OFF = ../../../utilities/utilities_mod.f90 ../../../common/types_mod.f90 ../../../preprocess/preprocess.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../time_manager/time_manager_mod.f90
clean: neat
	-rm -f .cppdefs $(OBJ) preprocess
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
preprocess: $(OBJ)
	$(LD) $(OBJ) -o preprocess  $(LDFLAGS)

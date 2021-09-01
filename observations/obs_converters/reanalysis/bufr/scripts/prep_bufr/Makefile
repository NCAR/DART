# Makefile created by mkmf $Id: mkmf 11289 2017-03-10 21:56:06Z hendric@ucar.edu $ 

include ../../../../../build_templates/mkmf.template


.DEFAULT:
	-touch $@
all: fixsys advance_time
fixsys: 
	@ (cd ../../../../../assimilation_code/modules/utilities; ./fixsystem $(FC) )
null_mpi_utilities_mod.o: ../../../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
time_manager_mod.o: ../../../../../assimilation_code/modules/utilities/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../../assimilation_code/modules/utilities/time_manager_mod.f90
parse_args_mod.o: ../../../../../assimilation_code/modules/utilities/parse_args_mod.f90 utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../../assimilation_code/modules/utilities/parse_args_mod.f90
types_mod.o: ../../../../../assimilation_code/modules/utilities/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../../assimilation_code/modules/utilities/types_mod.f90
utilities_mod.o: ../../../../../assimilation_code/modules/utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../../assimilation_code/modules/utilities/utilities_mod.f90
advance_time.o: ../../../../../assimilation_code/programs/advance_time/advance_time.f90 time_manager_mod.o utilities_mod.o parse_args_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../../../assimilation_code/programs/advance_time/advance_time.f90
./types_mod.f90: ../../../../../assimilation_code/modules/utilities/types_mod.f90
	cp ../../../../../assimilation_code/modules/utilities/types_mod.f90 .
./time_manager_mod.f90: ../../../../../assimilation_code/modules/utilities/time_manager_mod.f90
	cp ../../../../../assimilation_code/modules/utilities/time_manager_mod.f90 .
./advance_time.f90: ../../../../../assimilation_code/programs/advance_time/advance_time.f90
	cp ../../../../../assimilation_code/programs/advance_time/advance_time.f90 .
./utilities_mod.f90: ../../../../../assimilation_code/modules/utilities/utilities_mod.f90
	cp ../../../../../assimilation_code/modules/utilities/utilities_mod.f90 .
./parse_args_mod.f90: ../../../../../assimilation_code/modules/utilities/parse_args_mod.f90
	cp ../../../../../assimilation_code/modules/utilities/parse_args_mod.f90 .
./null_mpi_utilities_mod.f90: ../../../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
	cp ../../../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 .
SRC = ../../../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 ../../../../../assimilation_code/modules/utilities/time_manager_mod.f90 ../../../../../assimilation_code/modules/utilities/parse_args_mod.f90 ../../../../../assimilation_code/modules/utilities/types_mod.f90 ../../../../../assimilation_code/modules/utilities/utilities_mod.f90 ../../../../../assimilation_code/programs/advance_time/advance_time.f90
OBJ = null_mpi_utilities_mod.o time_manager_mod.o parse_args_mod.o types_mod.o utilities_mod.o advance_time.o
OFF = ../../../../../assimilation_code/modules/utilities/types_mod.f90 ../../../../../assimilation_code/modules/utilities/time_manager_mod.f90 ../../../../../assimilation_code/programs/advance_time/advance_time.f90 ../../../../../assimilation_code/modules/utilities/utilities_mod.f90 ../../../../../assimilation_code/modules/utilities/parse_args_mod.f90 ../../../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
clean: neat
	-rm -f .cppdefs $(OBJ) advance_time
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
advance_time: $(OBJ)
	$(LD) $(OBJ) -o advance_time  $(LDFLAGS)

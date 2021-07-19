# Makefile created by mkmf $Id$ 

include ../../../build_templates/mkmf.template


.DEFAULT:
	-touch $@
all: fixsys obs_seq_to_netcdf
fixsys: 
	@ (cd ../../../assimilation_code/modules/utilities; ./fixsystem $(FC) )
distributed_state_mod.o: ../../../assimilation_code/modules/utilities/distributed_state_mod.f90 null_mpi_utilities_mod.o types_mod.o ensemble_manager_mod.o null_win_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/distributed_state_mod.f90
location_io_mod.o: ../../../assimilation_code/location/utilities/location_io_mod.f90 types_mod.o netcdf_utilities_mod.o location_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/location/utilities/location_io_mod.f90
state_vector_io_mod.o: ../../../assimilation_code/modules/io/state_vector_io_mod.f90 adaptive_inflate_mod.o direct_netcdf_mod.o types_mod.o null_mpi_utilities_mod.o ensemble_manager_mod.o utilities_mod.o time_manager_mod.o io_filenames_mod.o model_mod.o state_structure_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/io/state_vector_io_mod.f90
ensemble_manager_mod.o: ../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o null_mpi_utilities_mod.o sort_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90
sort_mod.o: ../../../assimilation_code/modules/utilities/sort_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/sort_mod.f90
null_win_mod.o: ../../../assimilation_code/modules/utilities/null_win_mod.f90 types_mod.o ensemble_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/null_win_mod.f90
direct_netcdf_mod.o: ../../../assimilation_code/modules/io/direct_netcdf_mod.f90 types_mod.o options_mod.o ensemble_manager_mod.o time_manager_mod.o utilities_mod.o netcdf_utilities_mod.o null_mpi_utilities_mod.o state_structure_mod.o io_filenames_mod.o assim_model_mod.o model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/io/direct_netcdf_mod.f90
time_manager_mod.o: ../../../assimilation_code/modules/utilities/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/time_manager_mod.f90
netcdf_utilities_mod.o: ../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90
state_structure_mod.o: ../../../assimilation_code/modules/io/state_structure_mod.f90 utilities_mod.o obs_kind_mod.o types_mod.o netcdf_utilities_mod.o sort_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/io/state_structure_mod.f90
types_mod.o: ../../../assimilation_code/modules/utilities/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/types_mod.f90
assim_model_mod.o: ../../../assimilation_code/modules/assimilation/assim_model_mod.f90 types_mod.o time_manager_mod.o ensemble_manager_mod.o model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/assimilation/assim_model_mod.f90
obs_def_utilities_mod.o: ../../../observations/forward_operators/obs_def_utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/forward_operators/obs_def_utilities_mod.f90
schedule_mod.o: ../../../assimilation_code/modules/utilities/schedule_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/schedule_mod.f90
adaptive_inflate_mod.o: ../../../assimilation_code/modules/assimilation/adaptive_inflate_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o random_seq_mod.o ensemble_manager_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
default_model_mod.o: ../../../models/utilities/default_model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o netcdf_utilities_mod.o ensemble_manager_mod.o dart_time_io_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/utilities/default_model_mod.f90
io_filenames_mod.o: ../../../assimilation_code/modules/io/io_filenames_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o null_mpi_utilities_mod.o state_structure_mod.o ensemble_manager_mod.o netcdf_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/io/io_filenames_mod.f90
location_mod.o: ../../../assimilation_code/location/oned/location_mod.f90 types_mod.o utilities_mod.o random_seq_mod.o ensemble_manager_mod.o default_location_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/location/oned/location_mod.f90
random_seq_mod.o: ../../../assimilation_code/modules/utilities/random_seq_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/random_seq_mod.f90
obs_kind_mod.o: ../../../assimilation_code/modules/observations/obs_kind_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/observations/obs_kind_mod.f90
default_location_mod.o: ../../../assimilation_code/location/utilities/default_location_mod.f90 types_mod.o ensemble_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/location/utilities/default_location_mod.f90
options_mod.o: ../../../assimilation_code/modules/utilities/options_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/options_mod.f90
obs_sequence_mod.o: ../../../assimilation_code/modules/observations/obs_sequence_mod.f90 types_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/observations/obs_sequence_mod.f90
dart_time_io_mod.o: ../../../assimilation_code/modules/io/dart_time_io_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o netcdf_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/io/dart_time_io_mod.f90
model_mod.o: ../../../models/template/model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o netcdf_utilities_mod.o state_structure_mod.o ensemble_manager_mod.o dart_time_io_mod.o default_model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/template/model_mod.f90
obs_def_mod.o: ../../../observations/forward_operators/obs_def_mod.f90 types_mod.o utilities_mod.o location_mod.o obs_kind_mod.o assim_model_mod.o cov_cutoff_mod.o ensemble_manager_mod.o obs_def_utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/forward_operators/obs_def_mod.f90
cov_cutoff_mod.o: ../../../assimilation_code/modules/assimilation/cov_cutoff_mod.f90 types_mod.o utilities_mod.o location_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/assimilation/cov_cutoff_mod.f90
obs_seq_to_netcdf.o: ../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.f90 types_mod.o obs_sequence_mod.o obs_def_mod.o obs_kind_mod.o location_mod.o location_io_mod.o time_manager_mod.o schedule_mod.o utilities_mod.o netcdf_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.f90
null_mpi_utilities_mod.o: ../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
utilities_mod.o: ../../../assimilation_code/modules/utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assimilation_code/modules/utilities/utilities_mod.f90
./state_vector_io_mod.f90: ../../../assimilation_code/modules/io/state_vector_io_mod.f90
	cp ../../../assimilation_code/modules/io/state_vector_io_mod.f90 .
./default_model_mod.f90: ../../../models/utilities/default_model_mod.f90
	cp ../../../models/utilities/default_model_mod.f90 .
./state_structure_mod.f90: ../../../assimilation_code/modules/io/state_structure_mod.f90
	cp ../../../assimilation_code/modules/io/state_structure_mod.f90 .
./adaptive_inflate_mod.f90: ../../../assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
	cp ../../../assimilation_code/modules/assimilation/adaptive_inflate_mod.f90 .
./obs_kind_mod.f90: ../../../assimilation_code/modules/observations/obs_kind_mod.f90
	cp ../../../assimilation_code/modules/observations/obs_kind_mod.f90 .
./sort_mod.f90: ../../../assimilation_code/modules/utilities/sort_mod.f90
	cp ../../../assimilation_code/modules/utilities/sort_mod.f90 .
./netcdf_utilities_mod.f90: ../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90
	cp ../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90 .
./time_manager_mod.f90: ../../../assimilation_code/modules/utilities/time_manager_mod.f90
	cp ../../../assimilation_code/modules/utilities/time_manager_mod.f90 .
./random_seq_mod.f90: ../../../assimilation_code/modules/utilities/random_seq_mod.f90
	cp ../../../assimilation_code/modules/utilities/random_seq_mod.f90 .
./obs_seq_to_netcdf.f90: ../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.f90
	cp ../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.f90 .
./distributed_state_mod.f90: ../../../assimilation_code/modules/utilities/distributed_state_mod.f90
	cp ../../../assimilation_code/modules/utilities/distributed_state_mod.f90 .
./obs_def_utilities_mod.f90: ../../../observations/forward_operators/obs_def_utilities_mod.f90
	cp ../../../observations/forward_operators/obs_def_utilities_mod.f90 .
./schedule_mod.f90: ../../../assimilation_code/modules/utilities/schedule_mod.f90
	cp ../../../assimilation_code/modules/utilities/schedule_mod.f90 .
./utilities_mod.f90: ../../../assimilation_code/modules/utilities/utilities_mod.f90
	cp ../../../assimilation_code/modules/utilities/utilities_mod.f90 .
./default_location_mod.f90: ../../../assimilation_code/location/utilities/default_location_mod.f90
	cp ../../../assimilation_code/location/utilities/default_location_mod.f90 .
./obs_def_mod.f90: ../../../observations/forward_operators/obs_def_mod.f90
	cp ../../../observations/forward_operators/obs_def_mod.f90 .
./direct_netcdf_mod.f90: ../../../assimilation_code/modules/io/direct_netcdf_mod.f90
	cp ../../../assimilation_code/modules/io/direct_netcdf_mod.f90 .
./assim_model_mod.f90: ../../../assimilation_code/modules/assimilation/assim_model_mod.f90
	cp ../../../assimilation_code/modules/assimilation/assim_model_mod.f90 .
./ensemble_manager_mod.f90: ../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90
	cp ../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90 .
./obs_sequence_mod.f90: ../../../assimilation_code/modules/observations/obs_sequence_mod.f90
	cp ../../../assimilation_code/modules/observations/obs_sequence_mod.f90 .
./options_mod.f90: ../../../assimilation_code/modules/utilities/options_mod.f90
	cp ../../../assimilation_code/modules/utilities/options_mod.f90 .
./null_win_mod.f90: ../../../assimilation_code/modules/utilities/null_win_mod.f90
	cp ../../../assimilation_code/modules/utilities/null_win_mod.f90 .
./dart_time_io_mod.f90: ../../../assimilation_code/modules/io/dart_time_io_mod.f90
	cp ../../../assimilation_code/modules/io/dart_time_io_mod.f90 .
./cov_cutoff_mod.f90: ../../../assimilation_code/modules/assimilation/cov_cutoff_mod.f90
	cp ../../../assimilation_code/modules/assimilation/cov_cutoff_mod.f90 .
./model_mod.f90: ../../../models/template/model_mod.f90
	cp ../../../models/template/model_mod.f90 .
./null_mpi_utilities_mod.f90: ../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
	cp ../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 .
./types_mod.f90: ../../../assimilation_code/modules/utilities/types_mod.f90
	cp ../../../assimilation_code/modules/utilities/types_mod.f90 .
./io_filenames_mod.f90: ../../../assimilation_code/modules/io/io_filenames_mod.f90
	cp ../../../assimilation_code/modules/io/io_filenames_mod.f90 .
./location_io_mod.f90: ../../../assimilation_code/location/utilities/location_io_mod.f90
	cp ../../../assimilation_code/location/utilities/location_io_mod.f90 .
./location_mod.f90: ../../../assimilation_code/location/oned/location_mod.f90
	cp ../../../assimilation_code/location/oned/location_mod.f90 .
SRC = ../../../assimilation_code/modules/utilities/distributed_state_mod.f90 ../../../assimilation_code/location/utilities/location_io_mod.f90 ../../../assimilation_code/modules/io/state_vector_io_mod.f90 ../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90 ../../../assimilation_code/modules/utilities/sort_mod.f90 ../../../assimilation_code/modules/utilities/null_win_mod.f90 ../../../assimilation_code/modules/io/direct_netcdf_mod.f90 ../../../assimilation_code/modules/utilities/time_manager_mod.f90 ../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90 ../../../assimilation_code/modules/io/state_structure_mod.f90 ../../../assimilation_code/modules/utilities/types_mod.f90 ../../../assimilation_code/modules/assimilation/assim_model_mod.f90 ../../../observations/forward_operators/obs_def_utilities_mod.f90 ../../../assimilation_code/modules/utilities/schedule_mod.f90 ../../../assimilation_code/modules/assimilation/adaptive_inflate_mod.f90 ../../../models/utilities/default_model_mod.f90 ../../../assimilation_code/modules/io/io_filenames_mod.f90 ../../../assimilation_code/location/oned/location_mod.f90 ../../../assimilation_code/modules/utilities/random_seq_mod.f90 ../../../assimilation_code/modules/observations/obs_kind_mod.f90 ../../../assimilation_code/location/utilities/default_location_mod.f90 ../../../assimilation_code/modules/utilities/options_mod.f90 ../../../assimilation_code/modules/observations/obs_sequence_mod.f90 ../../../assimilation_code/modules/io/dart_time_io_mod.f90 ../../../models/template/model_mod.f90 ../../../observations/forward_operators/obs_def_mod.f90 ../../../assimilation_code/modules/assimilation/cov_cutoff_mod.f90 ../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.f90 ../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 ../../../assimilation_code/modules/utilities/utilities_mod.f90
OBJ = distributed_state_mod.o location_io_mod.o state_vector_io_mod.o ensemble_manager_mod.o sort_mod.o null_win_mod.o direct_netcdf_mod.o time_manager_mod.o netcdf_utilities_mod.o state_structure_mod.o types_mod.o assim_model_mod.o obs_def_utilities_mod.o schedule_mod.o adaptive_inflate_mod.o default_model_mod.o io_filenames_mod.o location_mod.o random_seq_mod.o obs_kind_mod.o default_location_mod.o options_mod.o obs_sequence_mod.o dart_time_io_mod.o model_mod.o obs_def_mod.o cov_cutoff_mod.o obs_seq_to_netcdf.o null_mpi_utilities_mod.o utilities_mod.o
OFF = ../../../assimilation_code/modules/io/state_vector_io_mod.f90 ../../../models/utilities/default_model_mod.f90 ../../../assimilation_code/modules/io/state_structure_mod.f90 ../../../assimilation_code/modules/assimilation/adaptive_inflate_mod.f90 ../../../assimilation_code/modules/observations/obs_kind_mod.f90 ../../../assimilation_code/modules/utilities/sort_mod.f90 ../../../assimilation_code/modules/utilities/netcdf_utilities_mod.f90 ../../../assimilation_code/modules/utilities/time_manager_mod.f90 ../../../assimilation_code/modules/utilities/random_seq_mod.f90 ../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.f90 ../../../assimilation_code/modules/utilities/distributed_state_mod.f90 ../../../observations/forward_operators/obs_def_utilities_mod.f90 ../../../assimilation_code/modules/utilities/schedule_mod.f90 ../../../assimilation_code/modules/utilities/utilities_mod.f90 ../../../assimilation_code/location/utilities/default_location_mod.f90 ../../../observations/forward_operators/obs_def_mod.f90 ../../../assimilation_code/modules/io/direct_netcdf_mod.f90 ../../../assimilation_code/modules/assimilation/assim_model_mod.f90 ../../../assimilation_code/modules/utilities/ensemble_manager_mod.f90 ../../../assimilation_code/modules/observations/obs_sequence_mod.f90 ../../../assimilation_code/modules/utilities/options_mod.f90 ../../../assimilation_code/modules/utilities/null_win_mod.f90 ../../../assimilation_code/modules/io/dart_time_io_mod.f90 ../../../assimilation_code/modules/assimilation/cov_cutoff_mod.f90 ../../../models/template/model_mod.f90 ../../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 ../../../assimilation_code/modules/utilities/types_mod.f90 ../../../assimilation_code/modules/io/io_filenames_mod.f90 ../../../assimilation_code/location/utilities/location_io_mod.f90 ../../../assimilation_code/location/oned/location_mod.f90
clean: neat
	-rm -f .cppdefs $(OBJ) obs_seq_to_netcdf
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
obs_seq_to_netcdf: $(OBJ)
	$(LD) $(OBJ) -o obs_seq_to_netcdf  $(LDFLAGS)

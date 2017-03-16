# Makefile created by mkmf $Id$ 

include ../../../mkmf/mkmf.template


.DEFAULT:
	-touch $@
all: fixsys gsi_to_dart
fixsys: 
	@ (cd ../../../mpi_utilities; ./fixsystem $(FC) )
mpisetup.o: ../../../observations/GSI2DART/mpisetup.f90 kinds.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/mpisetup.f90
assim_model_mod.o: ../../../assim_model/assim_model_mod.f90 types_mod.o location_mod.o time_manager_mod.o utilities_mod.o model_mod.o ensemble_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assim_model/assim_model_mod.f90
readsatobs.o: ../../../observations/GSI2DART/enkf/readsatobs.f90 kinds.o read_diag.o params.o radinfo.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/enkf/readsatobs.f90
params.o: ../../../observations/GSI2DART/params.f90 kinds.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/params.f90
state_structure_mod.o: ../../../io/state_structure_mod.f90 utilities_mod.o obs_kind_mod.o types_mod.o sort_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../io/state_structure_mod.f90
model_mod.o: ../../../models/wrf/model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o mpi_utilities_mod.o random_seq_mod.o obs_kind_mod.o ensemble_manager_mod.o sort_mod.o distributed_state_mod.o state_structure_mod.o module_map_utils.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/wrf/model_mod.f90
constants.o: ../../../observations/GSI2DART/enkf/constants.f90 kinds.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/enkf/constants.f90
obs_kind_mod.o: ../../../obs_kind/obs_kind_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_kind/obs_kind_mod.f90
no_cray_win_mod.o: ../../../distributed/no_cray_win_mod.f90 mpi_utilities_mod.o types_mod.o ensemble_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../distributed/no_cray_win_mod.f90
mpi_utilities_mod.o: ../../../mpi_utilities/mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../mpi_utilities/mpi_utilities_mod.f90
module_map_utils.o: ../../../models/wrf/module_map_utils.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/wrf/module_map_utils.f90
random_seq_mod.o: ../../../random_seq/random_seq_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../random_seq/random_seq_mod.f90
distributed_state_mod.o: ../../../distributed/distributed_state_mod.f90 mpi_utilities_mod.o types_mod.o ensemble_manager_mod.o no_cray_win_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../distributed/distributed_state_mod.f90
utilities_mod.o: ../../../utilities/utilities_mod.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../utilities/utilities_mod.f90
types_mod.o: ../../../common/types_mod.f90
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../common/types_mod.f90
dart_obs_seq_mod.o: ../../../observations/GSI2DART/dart_obs_seq_mod.f90 kinds.o params.o mpisetup.o obs_utilities_mod.o obs_sequence_mod.o obs_def_mod.o obs_kind_mod.o location_mod.o time_manager_mod.o radinfo.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/dart_obs_seq_mod.f90
mpi_readobs.o: ../../../observations/GSI2DART/enkf/mpi_readobs.f90 kinds.o radinfo.o readconvobs.o readsatobs.o readozobs.o mpisetup.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/enkf/mpi_readobs.f90
readconvobs.o: ../../../observations/GSI2DART/enkf/readconvobs.f90 kinds.o constants.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/enkf/readconvobs.f90
time_manager_mod.o: ../../../time_manager/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../time_manager/time_manager_mod.f90
kinds.o: ../../../observations/GSI2DART/enkf/kinds.f90
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/enkf/kinds.f90
radinfo.o: ../../../observations/GSI2DART/radinfo.f90 kinds.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/radinfo.f90
obs_utilities_mod.o: ../../../observations/utilities/obs_utilities_mod.f90 types_mod.o utilities_mod.o obs_def_mod.o obs_sequence_mod.o time_manager_mod.o location_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/utilities/obs_utilities_mod.f90
read_diag.o: ../../../observations/GSI2DART/enkf/read_diag.f90 kinds.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/enkf/read_diag.f90
obs_def_mod.o: ../../../obs_def/obs_def_mod.f90 types_mod.o utilities_mod.o location_mod.o assim_model_mod.o obs_kind_mod.o ensemble_manager_mod.o time_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_def/obs_def_mod.f90
gsi_to_dart.o: ../../../observations/GSI2DART/gsi_to_dart.f90 kinds.o mpisetup.o params.o radinfo.o mpi_readobs.o dart_obs_seq_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/gsi_to_dart.f90
readozobs.o: ../../../observations/GSI2DART/enkf/readozobs.f90 kinds.o params.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/GSI2DART/enkf/readozobs.f90
ensemble_manager_mod.o: ../../../ensemble_manager/ensemble_manager_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o mpi_utilities_mod.o sort_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../ensemble_manager/ensemble_manager_mod.f90
sort_mod.o: ../../../sort/sort_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../sort/sort_mod.f90
obs_sequence_mod.o: ../../../obs_sequence/obs_sequence_mod.f90 types_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_sequence/obs_sequence_mod.f90
location_mod.o: ../../../location/threed_sphere/location_mod.f90 types_mod.o utilities_mod.o random_seq_mod.o obs_kind_mod.o mpi_utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	../../../location/threed_sphere/location_mod.f90
./gsi_to_dart.f90: ../../../observations/GSI2DART/gsi_to_dart.f90
	cp ../../../observations/GSI2DART/gsi_to_dart.f90 .
./obs_sequence_mod.f90: ../../../obs_sequence/obs_sequence_mod.f90
	cp ../../../obs_sequence/obs_sequence_mod.f90 .
./mpisetup.f90: ../../../observations/GSI2DART/mpisetup.f90
	cp ../../../observations/GSI2DART/mpisetup.f90 .
./read_diag.f90: ../../../observations/GSI2DART/enkf/read_diag.f90
	cp ../../../observations/GSI2DART/enkf/read_diag.f90 .
./mpi_utilities_mod.f90: ../../../mpi_utilities/mpi_utilities_mod.f90
	cp ../../../mpi_utilities/mpi_utilities_mod.f90 .
./state_structure_mod.f90: ../../../io/state_structure_mod.f90
	cp ../../../io/state_structure_mod.f90 .
./location_mod.f90: ../../../location/threed_sphere/location_mod.f90
	cp ../../../location/threed_sphere/location_mod.f90 .
./readozobs.f90: ../../../observations/GSI2DART/enkf/readozobs.f90
	cp ../../../observations/GSI2DART/enkf/readozobs.f90 .
./readsatobs.f90: ../../../observations/GSI2DART/enkf/readsatobs.f90
	cp ../../../observations/GSI2DART/enkf/readsatobs.f90 .
./kinds.f90: ../../../observations/GSI2DART/enkf/kinds.f90
	cp ../../../observations/GSI2DART/enkf/kinds.f90 .
./ensemble_manager_mod.f90: ../../../ensemble_manager/ensemble_manager_mod.f90
	cp ../../../ensemble_manager/ensemble_manager_mod.f90 .
./obs_utilities_mod.f90: ../../../observations/utilities/obs_utilities_mod.f90
	cp ../../../observations/utilities/obs_utilities_mod.f90 .
./assim_model_mod.f90: ../../../assim_model/assim_model_mod.f90
	cp ../../../assim_model/assim_model_mod.f90 .
./sort_mod.f90: ../../../sort/sort_mod.f90
	cp ../../../sort/sort_mod.f90 .
./no_cray_win_mod.f90: ../../../distributed/no_cray_win_mod.f90
	cp ../../../distributed/no_cray_win_mod.f90 .
./distributed_state_mod.f90: ../../../distributed/distributed_state_mod.f90
	cp ../../../distributed/distributed_state_mod.f90 .
./constants.f90: ../../../observations/GSI2DART/enkf/constants.f90
	cp ../../../observations/GSI2DART/enkf/constants.f90 .
./module_map_utils.f90: ../../../models/wrf/module_map_utils.f90
	cp ../../../models/wrf/module_map_utils.f90 .
./time_manager_mod.f90: ../../../time_manager/time_manager_mod.f90
	cp ../../../time_manager/time_manager_mod.f90 .
./readconvobs.f90: ../../../observations/GSI2DART/enkf/readconvobs.f90
	cp ../../../observations/GSI2DART/enkf/readconvobs.f90 .
./utilities_mod.f90: ../../../utilities/utilities_mod.f90
	cp ../../../utilities/utilities_mod.f90 .
./obs_kind_mod.f90: ../../../obs_kind/obs_kind_mod.f90
	cp ../../../obs_kind/obs_kind_mod.f90 .
./params.f90: ../../../observations/GSI2DART/params.f90
	cp ../../../observations/GSI2DART/params.f90 .
./types_mod.f90: ../../../common/types_mod.f90
	cp ../../../common/types_mod.f90 .
./mpi_readobs.f90: ../../../observations/GSI2DART/enkf/mpi_readobs.f90
	cp ../../../observations/GSI2DART/enkf/mpi_readobs.f90 .
./radinfo.f90: ../../../observations/GSI2DART/radinfo.f90
	cp ../../../observations/GSI2DART/radinfo.f90 .
./model_mod.f90: ../../../models/wrf/model_mod.f90
	cp ../../../models/wrf/model_mod.f90 .
./obs_def_mod.f90: ../../../obs_def/obs_def_mod.f90
	cp ../../../obs_def/obs_def_mod.f90 .
./dart_obs_seq_mod.f90: ../../../observations/GSI2DART/dart_obs_seq_mod.f90
	cp ../../../observations/GSI2DART/dart_obs_seq_mod.f90 .
./random_seq_mod.f90: ../../../random_seq/random_seq_mod.f90
	cp ../../../random_seq/random_seq_mod.f90 .
SRC = ../../../observations/GSI2DART/mpisetup.f90 ../../../assim_model/assim_model_mod.f90 ../../../observations/GSI2DART/enkf/readsatobs.f90 ../../../observations/GSI2DART/params.f90 ../../../io/state_structure_mod.f90 ../../../models/wrf/model_mod.f90 ../../../observations/GSI2DART/enkf/constants.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../distributed/no_cray_win_mod.f90 ../../../mpi_utilities/mpi_utilities_mod.f90 ../../../models/wrf/module_map_utils.f90 ../../../random_seq/random_seq_mod.f90 ../../../distributed/distributed_state_mod.f90 ../../../utilities/utilities_mod.f90 ../../../common/types_mod.f90 ../../../observations/GSI2DART/dart_obs_seq_mod.f90 ../../../observations/GSI2DART/enkf/mpi_readobs.f90 ../../../observations/GSI2DART/enkf/readconvobs.f90 ../../../time_manager/time_manager_mod.f90 ../../../observations/GSI2DART/enkf/kinds.f90 ../../../observations/GSI2DART/radinfo.f90 ../../../observations/utilities/obs_utilities_mod.f90 ../../../observations/GSI2DART/enkf/read_diag.f90 ../../../obs_def/obs_def_mod.f90 ../../../observations/GSI2DART/gsi_to_dart.f90 ../../../observations/GSI2DART/enkf/readozobs.f90 ../../../ensemble_manager/ensemble_manager_mod.f90 ../../../sort/sort_mod.f90 ../../../obs_sequence/obs_sequence_mod.f90 ../../../location/threed_sphere/location_mod.f90
OBJ = mpisetup.o assim_model_mod.o readsatobs.o params.o state_structure_mod.o model_mod.o constants.o obs_kind_mod.o no_cray_win_mod.o mpi_utilities_mod.o module_map_utils.o random_seq_mod.o distributed_state_mod.o utilities_mod.o types_mod.o dart_obs_seq_mod.o mpi_readobs.o readconvobs.o time_manager_mod.o kinds.o radinfo.o obs_utilities_mod.o read_diag.o obs_def_mod.o gsi_to_dart.o readozobs.o ensemble_manager_mod.o sort_mod.o obs_sequence_mod.o location_mod.o
OFF = ../../../observations/GSI2DART/gsi_to_dart.f90 ../../../obs_sequence/obs_sequence_mod.f90 ../../../observations/GSI2DART/mpisetup.f90 ../../../observations/GSI2DART/enkf/read_diag.f90 ../../../mpi_utilities/mpi_utilities_mod.f90 ../../../io/state_structure_mod.f90 ../../../location/threed_sphere/location_mod.f90 ../../../observations/GSI2DART/enkf/readozobs.f90 ../../../observations/GSI2DART/enkf/readsatobs.f90 ../../../observations/GSI2DART/enkf/kinds.f90 ../../../ensemble_manager/ensemble_manager_mod.f90 ../../../observations/utilities/obs_utilities_mod.f90 ../../../assim_model/assim_model_mod.f90 ../../../sort/sort_mod.f90 ../../../distributed/no_cray_win_mod.f90 ../../../distributed/distributed_state_mod.f90 ../../../observations/GSI2DART/enkf/constants.f90 ../../../models/wrf/module_map_utils.f90 ../../../time_manager/time_manager_mod.f90 ../../../observations/GSI2DART/enkf/readconvobs.f90 ../../../utilities/utilities_mod.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../observations/GSI2DART/params.f90 ../../../common/types_mod.f90 ../../../observations/GSI2DART/enkf/mpi_readobs.f90 ../../../observations/GSI2DART/radinfo.f90 ../../../models/wrf/model_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../observations/GSI2DART/dart_obs_seq_mod.f90 ../../../random_seq/random_seq_mod.f90
clean: neat
	-rm -f .cppdefs $(OBJ) gsi_to_dart
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
gsi_to_dart: $(OBJ)
	$(MPILD) $(OBJ) -o gsi_to_dart  $(LDFLAGS)

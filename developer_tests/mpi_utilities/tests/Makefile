# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#

##
## Makefile for various mpi and netcdf test programs
##
## Usage: make <target>
##
# This Makefile directly includes the mkmf.template file used by the 
#  mkmf system.  If the mkmf.template is set wrong for this system, 
#  these builds will fail.  See ../../../build_templates directory
#  for various different compiler and platform files.

include ../../../build_templates/mkmf.template

# the include defines the following makefile variables:
# MPIFC =    # mpi fortran wrapper
# MPILD =    # mpi loader wrapper
# FC =       # fortran compiler
# LD =       # loader
# NETCDF =   # directory where include and lib for netcdf are located
# INCS =     # any additional include dirs needed
# FFLAGS += -I$(INCS)    # extra fortran compile flags
# LIBS =     # any additional libraries needed
# LDFLAGS += $(FFLAGS) -L$(NETCDF)/lib -lnetcdf  # loader flags

# if you want to test the C programs, you will need to
# define CC and MPICC here because they are not defined in
# the standard DART templates
CC    = gcc    # c compiler (e.g. gcc, icc, cc)
MPICC = mpicc  # mpi wrapper for cc

# one test uses the DART utilities, and several tests need
# to comment in or out blocks of code.

UTILS_DIR = ../../../assimilation_code/modules/utilities
FIXSYSTEM = $(UTILS_DIR)/fixsystem

# ----------------------------------------------------------------------------

SRC = ftest_f90.f90 ftest_mpi.f90 ftest_nc.f90 \
      ftest_onesided.f90 ftest_sendrecv.f90 \
      ctest.c ctest_mpi.c ctest_nc.c

ALLSRC = $(SRC) Makefile README runme_*

EXE = ftest_f90 ftest_mpi ftest_nc ftest_sendrecv ftest_onesided
ALLEXE = $(EXE) ctest ctest_mpi ctest_nc

FIXLIST = ftest_mpi.f90 

# ----------------------------------------------------------------------------
# default target
.DEFAULT_GOAL := help

## help: print help
.PHONY : help
help : 
	@sed -n 's/^##/   /p' $(MAKEFILE_LIST); echo ""

## all: build only what is necessary for dart
all:  $(EXE)

## everything: build all executables. Note set CC and MPICC to your C compiler
# if you are having problems, this might help diagnose
# what parts of the install are giving problems
everything: $(ALLEXE)

fixsys:
	$(FIXSYSTEM) $(FC) $(FIXLIST)


## check: run the basic executables interactively
check:  all
	./ftest_f90                   # this should run
	./ftest_nc                    # this should run
	mpirun -np 2 ./ftest_mpi      # this may not be allowed on a login node
	mpirun -np 2 ./ftest_onesided # this may not be allowed on a login node

## batch: submit the jobs to a batch queue
batch: ftest_mpi
	#bsub < runme           # comment this in for LSF
	#qsub runme             # comment this in for PBS
	./runme                 # anything else

# ----------------------------------------------------------------------------


# simple f90 test program without mpi or netcdf
ftest_f90: ftest_f90.f90
	$(FC) ftest_f90.f90 $(LDFLAGS) -o ftest_f90

# f90 test program that calls MPI functions
ftest_mpi: fixsys ftest_mpi.f90
	$(MPIFC) ftest_mpi.f90 $(LDFLAGS) -o ftest_mpi

# f90 test program that uses DART MPI utility routines
ftest_sendrecv: ftest_sendrecv.f90 types_mod.o utilities_mod.o time_manager_mod.o mpi_utilities_mod.o
	$(MPIFC) ftest_sendrecv.f90 types_mod.o utilities_mod.o time_manager_mod.o mpi_utilities_mod.o $(LDFLAGS) -o ftest_sendrecv

# f90 test program for large messages (greater than 2^31-1 elements)
ftest_sendrecv_big: ftest_sendrecv_big.f90 types_mod.o utilities_mod.o time_manager_mod.o mpi_utilities_mod.o
	$(MPIFC) ftest_sendrecv_big.f90 types_mod.o utilities_mod.o time_manager_mod.o mpi_utilities_mod.o $(LDFLAGS) -o ftest_sendrecv_big 

# test of MPI one-sided communication
ftest_onesided: ftest_onesided.f90 
	$(MPIFC) ftest_onesided.f90 $(LDFLAGS) -o ftest_onesided


types_mod.o: $(UTILS_DIR)/types_mod.f90
	$(MPIFC) $(FFLAGS) -c $<

utilities_mod.o: $(UTILS_DIR)/utilities_mod.f90
	$(MPIFC) $(FFLAGS) -c $<

time_manager_mod.o: $(UTILS_DIR)/time_manager_mod.f90
	$(MPIFC) $(FFLAGS) -c $<

mpi_utilities_mod.o: $(UTILS_DIR)/mpi_utilities_mod.f90
	cd $(UTILS_DIR); ./fixsystem $(FC)
	$(MPIFC) $(FFLAGS) -c $<

# f90 test program that calls netCDF functions
ftest_nc: ftest_nc.f90
	$(FC) ftest_nc.f90 $(LDFLAGS) -o ftest_nc



# DART uses no C code, but if you want to test whether your system has 
# working C interfaces for MPI or netCDF, you can use these programs.
# define CC and MPICC at the top of this Makefile, if needed, to reference your own C compiler

# simple c program without mpi or netcdf
ctest: ctest.c
	$(CC) ctest.c -o ctest

# c test program that calls MPI functions
ctest_mpi: ctest_mpi.c
	$(MPICC) ctest_mpi.c -o ctest_mpi

# c test program that calls netCDF functions
ctest_nc: ctest_nc.c
	$(CC) -I$(NETCDF)/include ctest_nc.c -L$(NETCDF)/lib -lnetcdf -o ctest_nc



# ----------------------------------------------------------------------------

## run_c: run test of C MPI
# not strictly a DART test; it has no C code.  But if there is a question
# of whether MPI was built without the F90 interface, this might help
# diagnose the problem.
run_c: ctest_mpi
	#bsub < runme_c      # comment this in for LSF
	#qsub runme_c        # comment this in for PBS
	./runme_c            # anything else


# ----------------------------------------------------------------------------
## sanity: checks f90 and mpi interactively and in batch mode
sanity:
	make clean
	make ftest_f90
	./ftest_f90
	@echo 'PASSED: able to compile an f90 program (no mpi)'
	make check
	@echo 'PASSED: able to compile and run an mpi program interactively'
	make batch
	@echo 'PASSED: able to compile and run an mpi program in batch'
	@echo 'PASSED: all sanity tests succeeded!'

# ----------------------------------------------------------------------------
## clean: remove autogenerated files
clean:
	rm -f *.o *.mod $(ALLEXE) out_* *testdata.nc *.orig

tar:
	tar -zcvf ~/comm.tar.gz $(ALLSRC)

untar:
	tar -zxvf ~/comm.tar.gz


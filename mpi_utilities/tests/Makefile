# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# Makefile for various mpi test and timing programs
#
# This Makefile directly 'include's the mkmf.template file used by the 
#  mkmf system.  If the mkmf.template is set wrong for this system, 
#  these builds will fail.  See ../../mkmf/mkmf.template.* for various
#  compiler and platform settings.

include ../../mkmf/mkmf.template

# the include defines the following makefile variables:
# MPIFC = 
# MPILD =
# FC =
# LD =
# NETCDF =
# INCS = 
# FFLAGS += -I$(INCS)  
# LIBS =
# LDFLAGS += -L$(NETCDF)/lib -lnetcdf

# ----------------------------------------------------------------------------

SRC = ftest_f90.f90 ftest_mpi.f90 ftest_nc.f90 \
      ftest_stop.f90 ftest_go.f90 ctest.c ctest_mpi.c ctest_nc.c

ALLSRC = $(SRC) Makefile README runme_*

EXE = ftest_f90 ftest_mpi ftest_nc
ALLEXE = $(EXE) ftest_stop ftest_go ctest ctest_mpi ctest_nc

# ----------------------------------------------------------------------------

# default target - build and test only what is necessary for dart
all:  $(EXE)


# all executables - if you are having problems, this might help diagnose
# what parts of the install are giving problems
everything: $(ALLEXE)


# run the basic executables interactively
check:  all
	./ftest_f90                   # this should run
	./ftest_nc                    # this should run
	mpirun -np 2 ./ftest_mpi      # this may not be allowed on a login node

# submit the jobs to a batch queue
batch: ftest_mpi
	#bsub < runme           # comment this in for LSF
	#qsub runme             # comment this in for PBS
	./runme                 # anything else

# check async=4 mode.  comment in the right line for your batch system.
async4: ftest_stop ftest_go
	#bsub < runme_async4      # comment this in for LSF
	#qsub runme_async4        # comment this in for PBS
	./runme_async4            # anything else

# ----------------------------------------------------------------------------


# simple f90 test program without mpi or netcdf
ftest_f90: ftest_f90.f90
	$(FC) ftest_f90.f90 $(LDFLAGS) -o ftest_f90

# f90 test program that calls MPI functions
ftest_mpi: ftest_mpi.f90
	$(MPIFC) ftest_mpi.f90 $(LDFLAGS) -o ftest_mpi

# f90 test program that calls netCDF functions
ftest_nc: ftest_nc.f90
	$(FC) ftest_nc.f90 $(LDFLAGS) -o ftest_nc


# test of named pipe (fifo) communication with mpi
ftest_stop: ftest_stop.f90
	$(MPIFC) ftest_stop.f90 $(LDFLAGS) -o ftest_stop

ftest_go: ftest_go.f90
	$(MPIFC) ftest_go.f90 $(LDFLAGS) -o ftest_go



# DART uses no C code, but if you want to test whether your system has 
# working C interfaces for MPI or netCDF, you can use these programs.
# Edit 'gcc' below, if needed, to reference your own C compiler.

# simple c program without mpi or netcdf
ctest: ctest.c
	gcc ctest.c -o ctest

# c test program that calls MPI functions
ctest_mpi: ctest_mpi.c
	mpicc ctest_mpi.c -o ctest_mpi

# c test program that calls netCDF functions
ctest_nc: ctest_nc.c
	gcc -I$(NETCDF)/include ctest_nc.c -L$(NETCDF)/lib -lnetcdf -o ctest_nc



# ----------------------------------------------------------------------------

# not strictly a DART test; it has no C code.  But if there is a question
# of whether MPI was built without the F90 interface, this might help
# diagnose the problem.
run_c: ctest_mpi
	#bsub < runme_c      # comment this in for LSF
	#qsub runme_c        # comment this in for PBS
	./runme_c            # anything else

# interactive test of the stop/go pair.   if you must run with a batch
# system, use the async4: target instead.
run_s:  ftest_stop ftest_go
	# needs batch script
	mpirun ./ftest_stop &
	sleep 5
	mpirun ./ftest_go


# ----------------------------------------------------------------------------

sanity:
	make ftest_f90
	./ftest_f90
	@echo 'PASSED: able to compile an f90 program (no mpi)'
	make check
	@echo 'PASSED: able to compile and run an mpi program interactively'
	make batch
	@echo 'PASSED: able to compile and run an mpi program in batch'
	@echo 'PASSED: all sanity tests succeeded!'

# ----------------------------------------------------------------------------

clean:
	rm -f *.o *.mod $(ALLEXE) out_* *testdata.nc

tar:
	tar -zcvf ~/comm.tar.gz $(ALLSRC)

untar:
	tar -zxvf ~/comm.tar.gz

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


# ----------------------------------------------------------------------------
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
#
# Makefile for various mpi test and timing programs
#
# This Makefile directly 'include's the mkmf.template file used by the 
#  mkmf system.  If the mkmf.template is set wrong for this system, 
#  these builds will fail.  See ../../mkmf/mkmf.template.* for various
#  compiler and platform settings.
#

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
	./ftest_f90
	./ftest_nc
	mpirun -np 2 ./ftest_mpi

# submit the jobs to a batch queue
batch: ftest_mpi
	#bsub < runme           # LSF
	#qsub runme             # PBS
	./runme                # everything else

# check async=4 mode.  comment in the right line for your batch system.
async4: ftest_stop ftest_go
	#bsub < runme_async4      # LSF
	#qsub runme_async4        # PBS
	./runme_async4           # everything else

# ----------------------------------------------------------------------------


# simple f90 test program without mpi or netcdf
ftest_f90: ftest_f90.f90
	$(FC) $(FFLAGS) ftest_f90.f90 $(LDFLAGS) -o ftest_f90

# f90 test program that calls MPI functions
ftest_mpi: ftest_mpi.f90
	$(MPIFC) $(FFLAGS) ftest_mpi.f90 $(LDFLAGS) -o ftest_mpi

# f90 test program that calls netCDF functions
ftest_nc: ftest_nc.f90
	$(FC) $(FFLAGS) ftest_nc.f90 $(LDFLAGS) -o ftest_nc


# test of named pipe (fifo) communication with mpi
ftest_stop: ftest_stop.f90
	$(MPIFC) $(FFLAGS) ftest_stop.f90 $(LDFLAGS) -o ftest_stop

ftest_go: ftest_go.f90
	$(MPIFC) $(FFLAGS) ftest_go.f90 $(LDFLAGS) -o ftest_go



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

# lam or mpich w/ or w/o batch system
#
# these rules are complicated by the fact that Makefiles run shell commands
# with sh syntax (because csh/tcsh cannot deal with multiline commands).
# but running bsub from sh inherits the sh environment and i have not
# been able to convince it to change to csh and read my .cshrc to get
# the modules and environment set up right before executing the program.
#
# now i am trying to make a csh shell script which invokes bsub -- but somehow
# it still is not working.  i do not like this approach because it because 
# adds yet another layer of indirection onto an already confusing setup, but
# it does work successfully.   this works fine from the command line, but not
# from inside the makefile.  
#
# batch systems are evil.
#

run_f: ftest_mpi
	if [ '$(BATCH)' = 'lsf' ] ; then \
	    ( /bin/tcsh -x -c < ./runme_lsf_f ) ;\
	    echo 'job submitted to LSF batch queue.' ;\
	elif [ '$(BATCH)' = 'pbs' ]; then \
	    qsub runme_pbs_f ;\
	    echo 'job submitted to PBS batch queue.' ;\
	elif [ '$(BATCH)' = 'nobatch' ]; then \
	    mpirun -v -np 2 -nolocal -machinefile ~/machinefile ./ftest_mpi > out_f 2>&1 ;\
	    echo 'done. output in file out_f' ;\
	else \
	    mpirun -np 2 ./ftest_mpi > out_f 2>&1 ;\
	    echo 'done. output in file out_f' ;\
	fi
	

run_s:  ftest_stop ftest_go
	# needs batch script
	mpirun ./ftest_stop &
	sleep 5
	mpirun ./ftest_go

run_c: ctest_mpi
	if [ '$(BATCH)' = 'lsf' ] ; then \
	    /bin/csh bsub < runme_lsf_c ;\
	    echo 'job submitted to LSF batch queue.' ;\
	elif [ '$(BATCH)' = 'pbs' ]; then \
	    qsub runme_pbs_c ;\
	    echo 'job submitted to PBS batch queue.' ;\
	elif [ '$(BATCH)' = 'nobatch' ]; then \
	    mpirun -v -np 2 -nolocal -machinefile ~/machinefile ./ctest > out_c 2>&1 ;\
	    echo 'done. output in file out_c' ;\
	else \
	    mpirun -np 2 ./ctest > out_c 2>&1 ;\
	    echo 'done. output in file out_c' ;\
	fi 
	

# ----------------------------------------------------------------------------

# these are specific to the lam mpi libs
setup:
	lamboot -ssi rpi tcp

takedown:
	lamhalt

status:
	lamtask


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



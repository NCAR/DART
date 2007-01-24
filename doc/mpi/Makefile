# ----------------------------------------------------------------------------
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# $Id$
#
# makefile for various mpi test and timing programs
#
#
# This Makefile has compile flags for a variety of systems.
# To select the right section, set either the Makefile variable
# or the shell environment variable DART_TARGET to one of the
# following values.
#
# In the Makefile, set it:
#   DART_TARGET = something_from_below
# 
# In your shell, set it:
#   (csh:)  setenv DART_TARGET something_from_below
#   (ksh:)  set DART_TARGET=something_from_below
#
#
# ----------------------------------------------------------------------------

# Strategy:  If the user knows what they are doing (or if there
# are multiple choices on a single machine), do not mess with the target
# if it is already set.  Otherwise, if we get the machine name and we know
# what works, set it here.  If we do not recognize the machine, error out.

ifeq ($(DART_TARGET),)

 # get the actual machine name and switch on that.
 HOST = $(shell uname -n)

 ifeq ($(HOST),dart.scd.ucar.edu)
   DART_TARGET=Linux_dart
 endif

 ifeq ($(HOST),fisher.scd.ucar.edu)
   DART_TARGET=Linux_dart
 endif

 ifeq ($(HOST),coral.ucar.edu)
   DART_TARGET=Linux_coral
 endif
 ifeq ($(HOST),cr0125en)
   DART_TARGET=Linux_coral
 endif
 ifeq ($(HOST),cr0122en)
   DART_TARGET=Linux_coral
 endif

 ifeq ($(HOST),bv1103en)
   DART_TARGET=Linux_bluevista
 endif
 ifeq ($(HOST),bv1203en)
   DART_TARGET=Linux_bluevista
 endif


 #ifeq ($(HOST),bangkok.cgd.ucar.edu)
 ifeq ($(HOST),bangkok)
   DART_TARGET=Linux_bangkok
 endif

 ifeq ($(HOST),ocotillo)
   DART_TARGET=Linux_ocotillo
 endif

 ifeq ($(HOST),ln0126en)
   DART_TARGET=Linux_lightning
 endif
 ifeq ($(HOST),ln0127en)
   DART_TARGET=Linux_lightning
 endif

 ifeq ($(HOST),cloud.scd.ucar.edu)
   DART_TARGET=Darwin_Absoft_Lam
 endif

 ifeq ($(HOST),moore.scd.ucar.edu)
   DART_TARGET=Darwin_Gfortran_Lam
 endif

 # If we still got here and the target it unset, it is a machine we do not
 # know.  The user is going to have to help us out.
 ifeq ($(DART_TARGET),)
   $(error "the environment or makefile variable DART_TARGET must be set")
  endif

endif


# if you want to know what we selected.  if this is annoying, comment it out.
$(warning "using $(DART_TARGET) for DART_TARGET")

# set default to be run on current node.  in the machine specific sections
# below, set which kind of batch queue system is on this machine.
BATCH = "interactive"

# mac laptop with absoft fortran (e.g. cloud)
ifeq ($(DART_TARGET),Darwin_Absoft_Lam)
MPI_HOME = /usr/local/lam
FC = f95
MPIFC = f95
# the only wrapper built by the 7.1.2 lam is g77 only.
# so we have to add the libs by hand.
#MPIFC = $(MPI_HOME)/bin/mpif77
#MPIFC = mpif90
INCS = -I$(MPI_HOME)/include 
LIBLIST = -llamf77mpi -lmpi -llam -lmpi
COMFLAGS = -X -bind_at_load
FMPIFLAGS = $(INCS)
FMPILIBS = -L$(MPI_HOME)/lib $(LIBLIST)
FMPILDFLAGS = $(COMFLAGS)
endif

# mac laptop with g95 open source (e.g. cloud)
ifeq ($(DART_TARGET),Darwin_G95_Openmpi)
MPI_HOME = /usr/local/openmpi
FC = g95
MPIFC = $(MPI_HOME)/bin/mpif90
endif

# mac laptop with gfortran open source (e.g. moore)
ifeq ($(DART_TARGET),Darwin_Gfortran_Lam)
MPI_HOME = /usr/local/lam
FC = gfortran
MPIFC = $(MPI_HOME)/bin/mpif90
endif

# ocotillo - Intel 8, mpich
#ifeq ($(DART_TARGET),Linux_Intel_8_Mpich)
ifeq ($(DART_TARGET),Linux_ocotillo)
# TODO: not tested yet.  no batch system; ssh only.
FC = ifort
MPIFC = mpif90
MPI_HOME = /usr/local/mpich
BATCH = nobatch
endif

# coral - Intel 9, mpich
#ifeq ($(DART_TARGET),Linux_Intel_9_Mpich)
ifeq ($(DART_TARGET),Linux_coral)
FC = ifort
MPIFC = mpif90
MPI_HOME = /coral/local/mpich64
BATCH = lsf
endif

# dart - Intel 9, lam mpi
#ifeq ($(DART_TARGET),Linux_Intel_9_Lam)
ifeq ($(DART_TARGET),Linux_dart)
FC = ifort
MPIFC = mpif77
MPI_HOME = /contrib/lam
endif

# bangkok - linux cluster, PGI compiler
#ifeq ($(DART_TARGET),Linux_PGI_9_Mpich)
ifeq ($(DART_TARGET),Linux_bangkok)
MPI_HOME = /usr/local/mpich-1.2.6-pgi-hpf-cc-5.2-4
MPIFC = $(MPI_HOME)/bin/mpif90
FC = pgf90
BATCH = pbs
endif

# lightning
#ifeq ($(DART_TARGET),Linux_PGI_9_Mpich)
ifeq ($(DART_TARGET),Linux_lightning)
MPI_HOME = /usr/local/mpich
MPIFC = $(MPI_HOME)/bin/mpif90
FC = pgf90
BATCH = lsf
endif

# bluevista
#ifeq ($(DART_TARGET),Linux_IBM)
ifeq ($(DART_TARGET),Linux_bluevista)
MPI_HOME = /usr/lib
MPIFC = mpxlf
FC = xlf
FFLAGS = -qsuffix=f=f90:cpp=F90 -q64 -qfree=f90
FMPIFLAGS = -qsuffix=f=f90:cpp=F90 -q64 -qfree=f90
BATCH = lsf
endif


# ----------------------------------------------------------------------------

# common to all platforms

OPT = -g
#OPT = -O2

SRC = ftest_f77.f ftest_f90.f90 ftest_mpi.f90 ctest.c ctest_mpi.c \
      driver.f90 commtest_mod.f90 ftest_nml.f90 

ALLSRC = $(SRC) Makefile README *.nml runme_*

EXE = ftest_f90 ftest_mpi commtest
ALLEXE = $(EXE) ftest_f77 ctest ctest_mpi ftest_nml \
         stest stest_mpi stest2a stest2b stest2a_mpi stest2b_mpi

# ----------------------------------------------------------------------------

# default target - build and test only what is necessary for dart
all:  $(EXE)


# all executables - if you are having problems, this might help diagnose
# what parts of the install are giving problems
everything: $(ALLEXE)


# run the basic executables interactively
check: 
	mpirun -np 2 ftest_mpi
	mpirun -np 4 commtest

# submit the jobs to a batch queue
batch:
	make run_f run_ct


# ----------------------------------------------------------------------------

# very simple f77 test program without mpi 
ftest_f77: ftest_f77.f
	$(FC) $(FFLAGS) ftest_f77.f $(FLIBS) -o ftest_f77

# very simple f90 test program without mpi 
ftest_f90: ftest_f90.f90
	$(FC) $(FFLAGS) ftest_f90.f90 -o $(FLIBS) ftest_f90

# f90 program which calls MPI functions
ftest_mpi: ftest_mpi.f90
	$(MPIFC) $(FMPIFLAGS) $(FMPILDFLAGS) ftest_mpi.f90 $(FMPILIBS) -o ftest_mpi

# very simple f90 test program with namelist input
ftest_nml: ftest_nml.f90
	$(FC) $(FFLAGS) ftest_nml.f90 -o $(FLIBS) ftest_nml

# test of named pipe communication, no mpi 
stest: stest.f90
	$(FC) $(FFLAGS) stest.f90 -o $(FLIBS) stest

# test of named pipe communication with mpi
stest_mpi: stest_mpi.f90
	$(MPIFC) $(FMPIFLAGS) stest_mpi.f90 -o $(FMPILIBS) stest_mpi

# test of paired programs, no mpi
stest2a: stest2a.f90
	$(FC) $(FFLAGS) stest2a.f90 -o $(FLIBS) stest2a

stest2b: stest2b.f90
	$(FC) $(FFLAGS) stest2b.f90 -o $(FLIBS) stest2b

# test of paired programs with mpi
stest2a_mpi: stest2a_mpi.f90
	$(MPIFC) $(FMPIFLAGS) stest2a_mpi.f90 -o $(FMPILIBS) stest2a_mpi

stest2b_mpi: stest2b_mpi.f90
	$(MPIFC) $(FMPIFLAGS) stest2b_mpi.f90 -o $(FMPILIBS) stest2b_mpi


# very simple c program without mpi
ctest: ctest.c
	gcc ctest.c -o ctest

# c program which calls MPI functions
ctest_mpi: ctest_mpi.c
	mpicc ctest.c -o ctest


# simple model for data assim communication pattern
commtest: commtest_mod.o driver.o
	$(MPIFC) $(FMPILDFLAGS) driver.o commtest_mod.o $(FMPILIBS) -o commtest

driver.o: driver.f90
	$(MPIFC) $(FMPIFLAGS) -c driver.f90

commtest_mod.o: commtest_mod.f90
	$(MPIFC) $(FMPIFLAGS) -c commtest_mod.f90



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

run_ct: commtest
	if [ '$(BATCH)' = 'lsf' ] ; then \
	    bsub < runme_lsf_ct ;\
	    echo 'job submitted to LSF batch queue.' ;\
	elif [ '$(BATCH)' = 'pbs' ]; then \
	    qsub runme_pbs_ct ;\
	    echo 'job submitted to PBS batch queue.' ;\
	elif [ '$(BATCH)' = 'nobatch' ]; then \
	    mpirun -v -np 4 -nolocal -machinefile ~/machinefile ./commtest > out_ct 2>&1 ;\
	    echo 'done. output in file out_ct' ;\
	else \
	    mpirun -np 4 ./commtest > out_ct 2>&1  ;\
	    echo 'done. output in file out_ct' ;\
	fi 
	
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
	

run_c: ctest_mpi
	if [ '$(BATCH)' = 'lsf' ] ; then \
	    bsub < runme_lsf_c ;\
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
	make ftest_f77
	./ftest_f77
	@echo 'PASSED: able to compile an f77 program (no mpi)'
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
	rm -f *.o *.mod $(ALLEXE) out_*

tar:
	tar -zcvf ~/comm.tar.gz $(ALLSRC)

untar:
	tar -zxvf ~/comm.tar.gz

ci:
	ci -l $(ALLSRC) 

co:
	co -l $(ALLSRC)

# see if there are diffs first; only get outputs on files which
# are different (this suppresses all the headers from identical files)
diff:
	@for i in $(ALLSRC) ;\
	do  \
	  ((rcsdiff -q $$i > /dev/null) || rcsdiff $$i) ;\
	done


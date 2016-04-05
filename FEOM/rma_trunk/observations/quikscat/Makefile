# This code is not protected by the DART copyright agreement.
# DART $Id$

# Makefile for Fortran QuikSCAT HDF readers
# Ted Lungu
#
# HDF libraries used to link the sample program
HDFDIR= /usr/local/hdf
HDFDIR= /contrib/hdf5-1.8.2_intel-10.1-64
HDFINC = $(HDFDIR)/include

LDFLAGS =  -L$(HDFDIR)/lib -lmfhdf -ldf -ljpeg -lz 
#******************************************************************
#       Uncomment next 4 lines if using SUN solaris machine
#******************************************************************
#MACHINE        =       SUN
#SYS_DEF        =       -lm -lnsl
#F77            =       f77
#FFLAGS         =       -w -e -cg92 -Nl99 -I$(HDFINC)

#******************************************************************
#    Uncomment next 4 lines if using Intel x86 Linux machine & g77
#******************************************************************
MACHINE        =       i586
SYS_DEF        =       -lm
F77            =       /contrib/hdf-4.2r4_gnu-4.1.2-64/bin/h4fc
FFLAGS         =       -ffixed-line-length-132 -O2 -I$(HDFINC)

#******************************************************************
#       Uncomment next 4 lines if using SUN OS
#******************************************************************
# MACHINE        =       SUN
# SYS_DEF        =       -lm
# F77            =       f77
# FFLAGS         =       -w -I$(HDFINC)

#******************************************************************
#       Uncomment next 3 lines if using SGI machine
#******************************************************************
# MACHINE        =       IRIX                                
# SYS_DEF        =       -lm
# F77            =       f77
# FFLAGS          =       -w -n32 -static -I$(HDFINC)

.KEEP_STATE:

read_qscat2b: read_qscat2b.o
	$(F77) read_qscat2b.o $(FFLAGS) $(LDFLAGS) $(SYS_DEF) -o $@

read_qscat2b.o: read_qscat2b.f
	$(F77) -c $(FFLAGS) read_qscat2b.f

clean:
	rm -f *.o

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


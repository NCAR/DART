SHELL=/bin/csh
#
FLNM= gswm
USER= hagan
USERUC= HAGAN
#CFT=	cft77
#CF=	f90
#FFLAGS= -Wf"-ez -o off"
#FFLAGS=
#
# Blackforest - http://www.scd.ucar.edu/docs/blackforest/index.html
#     compiling and building executables for single-thread codes
FC      = xlf_r
#FFLAGS  = -qarch=pwr3 -qrealsize=8 -O3 -qstrict
# as suggested by Dick Valent; email sent Dec 27, 2001:
FFLAGS  = -q32 -qrealsize=8 -g \
	  -qalias=noaryovrlp -O3 -qstrict \
          -qflttrap=ov:zero:invalid:enable \
	  -qsigtrap=xl__trcedump

#LDFLAGS=  -L /usr/local/lib32/r8i4/liblapack.a
# as suggested by Dick Valent; email sent Dec 27, 2001:
LDFLAGS=  -L /usr/local/lib32/r8i4 -lfitpack -llinpack \
          -bloadmap:loadmap -bmaxdata:0X70000000
#
# Blackforest:
COMMAND=$(FLNM)
PROFLIB=
MSPATH=/$(USERUC)/$(FLNM)
LIST = $(FLNM).list
#
SRCS=\
abcr1_2.f abcrn.f ali.f alpbet2.f amphz.f atmos5.f atsm.f \
setatmos.f msise90.f conv.f cool.f eddy_gs.f htvec.f det3.f \
forit.f geogm.f gmaddm.f gmaddv.f heat92.f \
ionden2.f ionos3.f latvecs.f seteddy_gs.f \
main_2.f mprdd.f mprddv.f polar.f prnt.f \
psi.f qsf.f semian.f sol.f tvarf2.f tvef1.f \
varx.f w.f yonii.f zeta.f driver.f setheat.f fric.f \
seth2o.f \
hwm93.f h_diff.f
#
OBJS=\
abcr1_2.o abcrn.o ali.o alpbet2.o amphz.o atmos5.o atsm.o \
setatmos.o msise90.o conv.o cool.o eddy_gs.o htvec.o det3.o \
forit.o geogm.o gmaddm.o gmaddv.o heat92.o \
ionden2.o ionos3.o latvecs.o seteddy_gs.o \
main_2.o mprdd.o mprddv.o polar.o prnt.o \
psi.o qsf.o semian.o sol.o tvarf2.o tvef1.o \
varx.o w.o yonii.o zeta.o driver.o setheat.o fric.o \
seth2o.o \
hwm93.o h_diff.o
#
INCH = 
INCF = 
#
#.KEEP_STATE:
#
# Make executable:
#

$(COMMAND):	$(OBJS)
	$(FC) $(FFLAGS) -o $(COMMAND) $(OBJS) $(LDFLAGS)
#
#
# Make profiled version of shavano executable:
# WARNING:  DO NOT PROFILE AND FLOWTRACE AT THE SAME TIME
#
#$(COMMAND).prof:	$(OBJS)
#	$(LDR) $(LDFLAGS) -o $(COMMAND) $(OBJS) $(PROFLIB) $(LIBS)
#
#
# flint is fortran source code analyzer: (jkh)
#
#lint:
#	flint -afgsmtx -Ttrim -Xnoline,tabular -S$(FLNM) $(SRCS)
#clean:
#	rm M* *.f *.h *.c *.tar *.inp *.out *.o *.ymp *.j90 *.jse gmeta core


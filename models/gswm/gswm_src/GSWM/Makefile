SHELL=/bin/csh
#
FLNM= gswm
USER= hackney
USERUC= HACKNEY
#CFT=	cft77
CF=	f90
#FFLAGS= -Wf"-ez -o off"
FFLAGS=
#LDR=	segldr
#LDFLAGS= -D dupentry=ignore -i /u2/$(USER)/segdir
LDFLAGS=  -L /usr/lib,/usr/local/lib,/lib \
	-l mv1,ncarg,ncarg_gks,ncarg_c,ncarg_loc,X11,ncarm, \
	-l ncaro,mss,alfpack,ecmfft,fitpack,fishpak,fftpack,hpf
#
# Shavano:
#COMMAND=$(FLNM).ymp
# Paiute:
#COMMAND=$(FLNM).j90
# Ouray:
COMMAND=$(FLNM).jse
PROFLIB=
MSPATH=/$(USERUC)/$(FLNM)
LIST = $(FLNM).list
#
SRCS=\
abcr1_2.f abcrn.f ali.f alpbet2.f amphz.f atmos5.f atsm.f \
setatmos.f msise90.f conv.f cool.f eddy_gs.f htvec.f det3.f \
forit.f geogm.f getlate.f gmaddm.f gmaddv.f heat92.f \
honglatent.f ionden2.f ionos3.f latvecs.f seteddy_gs.f \
main_2.f mprdd.f mprddv.f northv.f polar.f prnt.f geo.f \
psi.f qsf.f semian.f sol.f tvarf2.f tvef1.f wind2.f \
varx.f w.f yonii.f zeta.f driver.f setheat.f fric.f \
hoff.f hoffi.f seth2o.f latent.f latent_nm.f hme11.f eddy.f \
amphzgeo.f h2onm.f setgci.f seth2onm.f set_euv.f euvheat.f \
heat.f h_diff.f
#
OBJS=\
abcr1_2.o abcrn.o ali.o alpbet2.o amphz.o atmos5.o atsm.o \
setatmos.o msise90.o conv.o cool.o eddy_gs.o htvec.o det3.o \
forit.o geogm.o getlate.o gmaddm.o gmaddv.o heat92.o \
honglatent.o ionden2.o ionos3.o latvecs.o seteddy_gs.o \
main_2.o mprdd.o mprddv.o northv.o polar.o prnt.o geo.o \
psi.o qsf.o semian.o sol.o tvarf2.o tvef1.o wind2.o \
varx.o w.o yonii.o zeta.o driver.o setheat.o fric.o \
hoff.o hoffi.o seth2o.o latent.o latent_nm.o hme11.o eddy.o \
amphzgeo.o h2onm.o setgci.o seth2onm.o set_euv.o euvheat.o \
heat.o h_diff.o
#
INCH = 
INCF = 
#
.KEEP_STATE:
#
# Make executable: uncomment to run (jkh)
#
$(COMMAND):	$(OBJS)
	$(CF) -o $(COMMAND) $(OBJS) $(LDFLAGS)
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



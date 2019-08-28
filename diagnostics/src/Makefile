
include ../../../../build_templates/mkmf.template

CPP=cpp -P -C

######################################################
# Preprocessor definition

CPP_DEFS= -DPETSC -Duse_tracer_fct -Duse_ice -Duse_ice_fct -Dallow_diag -Dallow_calcmeans -Duse_fullfreesurf

# modules
MODULES=oce_modules.o            \
        utilities.o              \
        gen_modules_config.o     \
        random_perturbation.o    \
        gen_modules_clock.o      \
        fesom_forcing_diag.o     \
        fesom_obs_diag.o         \
        fesom_dart_diag.o        \
        fesom_ocean_diag.o

# objects
OBJECTS=oce_dens_press.o  \
        gen_input.o       \
        oce_mesh_setup.o  \
        mesh_read.o       \
        fesom_post_main.o

######################################################

EXE = ../fesom_post_mesh.x

$(EXE) : $(MODULES) $(OBJECTS)
	rm -f $@
	$(FC)  $(FFLAGS) -o $@ $(MODULES) $(OBJECTS) $(LIBS) #-lessl

.SUFFIXES: .F90 .o

.F90.o :
	$(CPP) $(CPP_DEFS) $(INCS) < $*.F90 > $*.f
	$(FC) $(FFLAGS) -c $*.f 
	rm -f  $*.f

clean :
	rm -f *.o *.mod *~ *.f $(EXE)

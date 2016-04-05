#
# There are three files that need to be present in this directory if you want to
# use this Makefile instead of using  ../work/mkmf_driver
#
# ln -s ../../../types_mod/types_mod.f90 types_mod.F
# ln -s ../../../utilities/utilities_mod.f90 utilities_mod.F
# ln -s ../../../time_manager/time_manager_mod.f90 time_manager_mod.F
#



NETCDF = /contrib/netcdf-3.6.0p1-pgcc-pgf90.6.0-5
FC = pgf90
FCFLAGSFREE = -g -C -byteswapio -Mfree -r8 -Ktrap=fp -pc 64
#LINKFLAGS = -g
#FC = ifc
#FCFLAGSFREE = -FR -g -C -Vaxlib -r8
#LINKFLAGS = -g -Vaxlib
LIB= -L$(NETCDF)/lib -lnetcdf
INCL= $(NETCDF)/include 
RM = /bin/rm -f

.SUFFIXES: .F .f90 .o 

.F.o:
	$(FC) -c $(FCFLAGSFREE) -I$(INCL) $*.F

.f90.o:
	$(FC) -c $(FCFLAGSFREE) -I$(INCL) $*.f90


MODULES =  module_wrf_error.o module_namelist.o   module_ideal.o \
module_initialize.o  \
module_model_constants.o module_soil_pre.o \
module_bl_ysu.o module_bl_mrf.o  module_bl_myjpbl.o \
module_sf_sfclay.o module_sf_slab.o  module_sf_force_restore.o \
module_sf_noahlsm.o  \
module_sf_sfcdiags.o module_sf_myjsfc.o module_sf_ruclsm.o \
module_pbl_driver.o  module_surface_driver.o  \
module_nr_procedures.o module_luse_init.o \
time_manager_mod.o types_mod.o utilities_mod.o \
module_snd_init_and_bc.o module_wrf_init_and_bc.o \
module_getsm.o module_init_soil_ideal.o module_init_soil_real.o \
module_interpolations.o	module_wrf.o 


OBJS = driver.o tridiag.o

driver.x : $(MODULES) $(OBJS)
	$(FC) $(LINKFLAGS) -I$(INCL) $(MODULES) $(OBJS) $(LIB) -o $@

driver.o: module_wrf.o
module_wrf.o:  module_model_constants.o module_namelist.o \
module_ideal.o module_initialize.o module_soil_pre.o \
module_snd_init_and_bc.o module_wrf_init_and_bc.o \
module_init_soil_ideal.o module_init_soil_real.o \
module_surface_driver.o module_pbl_driver.o module_sf_myjsfc.o \
module_bl_myjpbl.o module_sf_sfclay.o module_bl_mrf.o \
module_bl_ysu.o module_sf_noahlsm.o module_sf_ruclsm.o \
module_sf_force_restore.o \
module_luse_init.o module_getsm.o
module_surface_driver.o: module_model_constants.o module_sf_sfclay.o \
module_sf_myjsfc.o module_sf_noahlsm.o module_sf_ruclsm.o \
module_sf_sfcdiags.o module_sf_slab.o module_ideal.o \
module_sf_force_restore.o
module_initialize.o: module_namelist.o   module_interpolations.o
module_namelist.o: module_model_constants.o
module_wrf_init_and_bc.o: time_manager_mod.o module_model_constants.o \
module_namelist.o
module_snd_init_and_bc.o: time_manager_mod.o module_model_constants.o \
module_namelist.o
time_manager_mod.o: types_mod.o utilities_mod.o

clean:
	$(RM) *.o *.x *.mod 




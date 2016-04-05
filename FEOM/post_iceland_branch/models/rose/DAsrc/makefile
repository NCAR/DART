GOAL = DArose 
OBJ = params.mod.o \
      rose.mod.o \
      dynam.mod.o \
      chem0d.mod.o \
      chem.mod.o \
      phys.mod.o \
      airglow.mod.o \
      ncdf.mod.o \
      trop.mod.o \
      output.mod.o \
      diagnostic.mod.o \
      a3d.chem.o \
      chemheat.o \
      chembc_intrp.o \
      chem_gs.o \
      chem_nr.o \
      cool_ir.z.o \
      density.o \
      dmolec.o \
      dynamlbc.o \
      fft.o \
      filter.o \
      gravity.o \
      hadjust_fix.o \
      heat_uv.o \
      hines_driver.o \
      hines.o \
      hydrost.o \
      moldiff.o \
      msetfix.o \
      msetvar.o \
      msis90.o \
      nmc_trop.o \
      nr_solve.o \
      phot.o \
      radiation.o \
      rates.o \
      select0d.o \
      shapirx.o \
      shapiry.o \
      shapir_w.o \
      svdiffbc.o \
      svdiffch.o \
      tendency.o \
      tidalinp.o \
      tinterp_dyn.o \
      tinterp_rad.o \
      total.o \
      transport.o \
      update0d.o \
      vertdiff_c.o \
      vertdiff_t.o \
      wcont.o \
      zenith.o

# compiler options for model runs

# Portran Group Compiler Settings
# FFLAGS= -c -Ktrap=divz,ovf -fast -I/usr/include -Minfo
  FC= pgf90
  FFLAGS= -c -Ktrap=divz,ovf -fast -I/usr/local/netcdf-pgi-hpf-cc/include 
  LIBS= -L/usr/local/netcdf-pgi-hpf-cc/lib -lnetcdf 

# Intel Compiler Settings
#FC= ifc  # Portland compiler
#FFLAGS= -c -pc64 -zero -O0 -I/opt/Intel/include 
#LIBS= -L/opt/Intel/lib -lnetcdf

# dont know which compiler uses these flags 
# FFLAGS= -C -c -Ktrap=divz,ovf  # <<== Very slow - for degugging only
# flags for timing
# FFLAGS= -c -fast -Mvect=prefetch -Ktrap=divz,ovf -Mprof=func -Minfo
# LFLAGS= -Mprof=func
LFLAGS = 

.SUFFIXES : .f .o

.f.o :
	$(FC) $(FFLAGS) $<

#  link -- with netcdf libraries
program : $(OBJ)
	$(FC) $(LFLAGS) -o $(GOAL) $(OBJ) $(LIBS) 

clean:
	rm -v *.o *.mod $(GOAL)


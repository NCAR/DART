#!/bin/csh
#
# biogeochem/CNBalanceCheckMod.F90
# biogeophys/SurfaceRadiationMod.F90
# biogeophys/PhotosynthesisMod.F90
# biogeophys/CanopyFluxesMod.F90
# cpl/lnd_import_export.F90

foreach FILE ( \
               biogeochem/CNBalanceCheckMod.F90 \
               biogeophys/SurfaceRadiationMod.F90 \
               biogeophys/PhotosynthesisMod.F90 \
               biogeophys/CanopyFluxesMod.F90 )

   set OLD = /glade/work/thoar/CESM/cesm2.1.0/components/clm/src/$FILE
   set NEW = /glade/work/thoar/CESM/my_cesm_sandbox/components/clm/src/$FILE

#  diffuse $OLD $FILE
#  diffuse $OLD $FILE $NEW
   diffuse      $FILE $NEW

end

# The lnd_import_export.F90 is not in the same directory as it was ... I am not sure it is needed

   set FILE = cpl/lnd_import_export.F90
   set OLD = /glade/work/thoar/CESM/cesm2.1.0/components/clm/src/cpl/lnd_import_export.F90
   set NEW = /glade/work/thoar/CESM/my_cesm_sandbox/components/clm/src/cpl/mct/lnd_import_export.F90

   diffuse $OLD $FILE
   diffuse $OLD $FILE $NEW
   diffuse      $FILE $NEW

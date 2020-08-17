#!/bin/csh -f
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# ---------------------------------------------------------------------
# Converts 'old' CLM restart files to whatever resolution you like.
# I ran a default case for the time period and compset of interest and
# then used one of the CLM restart files as the 'template' below.
#
# The source CLM restart files were from some previous run.
# ---------------------------------------------------------------------
#
#BSUB -J interpinic
#BSUB -o interpinic.%J.log
#BSUB -P P8685nnnn
#BSUB -q geyser
#BSUB -N -u ${USER}@ucar.edu
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -W 3:00

set EXEDIR = ~thoar/cesm1_1_1/clm/tools/interpinic
set SRCDIR = /glade/scratch/thoar/DART_POP_RESTARTS/2004-01-01-00000
set OUTDIR = /glade/scratch/thoar/DART_POP_RESTARTS/

cd $OUTDIR

@ instance = 1

while ($instance <= 30)

   echo `date`" converting CLM restart file for ensemble member $instance ... "

   set file1 = `printf b40.20th.005_ens%02d.clm2.r.2004-01-01-00000.nc $instance`
   set file2 = `printf cesm_test.clm2_%04d.r.2004-01-01-00000.nc $instance`

   \cp -f clm_template_restart_file.nc $file2

   ${EXEDIR}/interpinic -i ${SRCDIR}/$file1 -o $file2

   @ instance ++
end

exit 0



#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download


# Compiles pert_sounding.f90 for ocotillo

#ifort -c pert_sounding_module.f90
#ifort pert_sounding.f90 -o pert_sounding pert_sounding_module.o
#\rm ./pert_sounding_mod.mod



# Compiles pert_sounding.f90 for bluefire

xlf -c pert_sounding_module.f90
xlf pert_sounding.f90 -o pert_sounding pert_sounding_module.o
rm ./pert_sounding_mod.mod

exit $status



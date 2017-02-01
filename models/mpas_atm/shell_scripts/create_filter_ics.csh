#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
##############################################################################################
#
# create_filter_ics.csh
# To create initial ensemble in the DART binary format,
# run model_to_dart over the mpas_init.nc for each member at the initial cycle time.
# 
# So-Young Ha (MMM/NCAR)
#
##############################################################################################
# USER SPECIFIED PARAMETERS
##############################################################################################
set  expname = MPAS_DART_test
set init_cyc = 2008080112						# initial cycle time
set   f_mpas = mpas_init.2008-08-01_12.00.00.nc				# mpas forecast file name
set dir_fcst = /glade/scratch/syha/MPAS/INIT_PERT			# ensemble forecast directory
set  dir_out = /glade/scratch/syha/MPAS_DART/$expname/${init_cyc}	# output directory
set dir_dart = /glade/scratch/syha/DART/branch_dev/models/mpas_atm/work	# dart executables

set n = 96		# ensemble size
set frst = filter_ics	# restart_in_file_name in &filter_nml in input.nml 
##############################################################################################
# END OF USER SPECIFIED PARAMETERS
##############################################################################################

set  infn = `grep   model_analysis_filename input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`
set outfn = `grep model_to_dart_output_file input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`

foreach f ( $infn $frst model_to_dart.log )
  if(-e $f) \rm -f $f
end
if(! -d $dir_out) mkdir -p $dir_out
cd $dir_out

if(! -e input.nml) \cp -f $dir_dart/input.nml .

set i = 1
while ( $i <= $n )
  set fsrc = $dir_fcst/ENS_$i/$f_mpas
  if(! -e $fsrc) then
     echo $fsrc does not exist. Stop.
     exit
  endif
  echo ln -sf $fsrc $infn
       ln -sf $fsrc $infn
  $dir_dart/model_to_dart >>&! model_to_dart.e$i.log

  if(! -e $outfn) then
     Error in model_to_dart for ensemble $i. Stop.
     exit
  endif
  set icnum = `echo $i + 10000 | bc | cut -b2-5`
  cat $outfn >! $dir_out/$frst.${icnum}

  set ff = `echo $infn | cut -d . -f1`
  ln -sf $fsrc ${ff}.e${i}.nc
@ i++
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


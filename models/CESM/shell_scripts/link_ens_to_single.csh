#!/bin/tcsh

# Script to make an ensemble of links to a single CESM restart set
# for use by the first cycle of a CESM+DART assimilation.
# Run this in the directory where the single restart set lives,
# or change the links below to have the right directory name(s).

# Set the $case part of the filename to which you want to link.
set refcase =  CAM6_spinup_sst.25
# Set the date and time of the filename to which you want to link.
set cesm_time = 2013-06-15-00000
# Set the ensemble size
@ num_instances = 80

@ inst=1
while ($inst <= $num_instances)

   set inst_string = `printf _%04d $inst`

   # echo "Staging initial files for instance $inst of $num_instances"

   ln -s ./${refcase}.clm2.r.${cesm_time}.nc ./${refcase}.clm2${inst_string}.r.${cesm_time}.nc 
   ln -s ./${refcase}.cice.r.${cesm_time}.nc ./${refcase}.cice${inst_string}.r.${cesm_time}.nc 
   ln -s ./${refcase}.cam.i.${cesm_time}.nc  ./${refcase}.cam${inst_string}.i.${cesm_time}.nc  

   # If you are using a river runoff model, you must specify an initial file here
   # ln -s ./${refcase}.rtm.r.${cesm_time}.nc ./${refcase}.rtm${inst_string}.r.${cesm_time}.nc
   # OR
   ln -s ./${refcase}.mosart.r.${cesm_time}.nc ./${refcase}.mosart${inst_string}.r.${cesm_time}.nc

   # If CISM is active
   ln -s ./${refcase}.cism.r.${cesm_time}.nc ./${refcase}.cism${inst_string}.r.${cesm_time}.nc
   
  @ inst ++
end

exit
   

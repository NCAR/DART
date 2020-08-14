#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Extract the mean copy for each field from each of the timeslots on
# a series of files (Posterior usually) and put each timeslot into 
# CAM initial file format.
# What about the CLM initial files?  Can't average some fields on them.
#    Does it make a difference to short forecasts?

# set echo verbose

# P[oste]rior_core.nc files will have copies 1,2,83,84, but ordered as 83,84,1,2
# ncks will call copy 83 'copy 1', so those files will require 
# dim_val = 3

set ms_file   = $1
set local_dir = $2
set kind      = $3
set dim       = $4
set element1  = $5
set yrmoday   = $6

set num_times  = 4
set hours      = (06 12 18 24)
set ens_member = 1

if (! -d $local_dir) mkdir $local_dir
cd $local_dir

if (! -e ${kind}.nc) then
   if (! -e ${kind}_Diag.nc) then
      if (! -e diagnostics.tar) then
         if (! -e diagnostics.tar.gz) then
            echo "starting msrcp at `date`"
            msrcp mss:$ms_file . 
         endif
         echo "starting gunzip at `date`"
# ? skip this and add -z to tar -x below?  That would short-circuit some of
#   the failure recovery I've built in, but in this context I can't use it.
         gunzip diagnostics.tar.gz
      endif
      echo "starting tar -x at `date`"
      tar -x -f diagnostics.tar ${kind}_Diag.nc
      if ($status == 0) rm diagnostics.tar &
   endif

   set time = 1
   while ($time <= $num_times)
      set out_name = cam_init_${yrmoday}$hours[$time].nc
      # Template into which we'll stuff the mean fields
#      cp ../caminput.nc $out_name
#     This could be a mv, if it's reliable enough
      cp cam_init_memb${ens_member}_H$hours[$time].nc    $out_name

      # ncks makes the time/record dimension of the out_file to be the largest 
      # of the sizes of the in_ and out_file, if out_file pre-exists.  
      # Prevent this by extracting one timeslot to a new file.
      # (This had -x -v copy in my first successful test, but that's redundant(?) 
      # with the ncwa, where'd I'd like to handle all the 'copy' purging.
      ncks -F -A -d ${dim},${element1} -d time,$time ${kind}_Diag.nc time_copy_slab.nc 
      if ($status != 0) then
         echo "timeslot $time not found on ${kind}_Diag.nc"
         exit
      endif
      # Exclude the copy dimension and variable from new file, 
      # along with mis-dimensioned P0
      # ncwa will change the rank; use it to remove dimension 'copy'
      ncwa -x -v P0,CopyMetaData -a ${dim} time_copy_slab.nc avgd_copy_out.nc
      # Permute the remaining dimensions into CAM; (time, lev, [s]lat, [s]lon)
      #                                  from DART (time, lat, [s]lon, [s]lev)
      # Requires 3 calls to handle the 3 (staggered) grids
      ncpdq -a lev,lat,lon  avgd_copy_out.nc re-orderedT.nc
      ncpdq -a lev,slat,lon re-orderedT.nc   re-orderedUS.nc
      ncpdq -a lev,lat,slon re-orderedUS.nc  re-orderedVS.nc
      # Finally, stuff the contents into the CAM initial file template
      # -A causes variables to be overwritten if there is a name "conflict",
      # without querying the user/script.
      # -x -v $dim removes the *variable* called copy.  (On blueice) it didn't work
      # to try to remove it in the ncwa command.
      ncks -A -x -v ${dim} re-orderedVS.nc $out_name

      @ time++
      rm time_copy_slab.nc avgd_copy_out.nc re-order*
   end
   cd ..

endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


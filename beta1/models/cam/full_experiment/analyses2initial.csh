#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Extract the mean copy for each state vector field from each of the timeslots 
# on a file (Posterior usually) and put each timeslot into CAM initial file format.
# The non-state fields are ensemble averages of all the caminput_#.nc from the same 
# timeslot.

# Unlike previous versions of this script, CLM ensemble files are not averaged
# because it is unclear whether it can ever be made to work.

# set echo verbose


# Called from auto_diag2ms_LSF.csh with
#   ../../analyses2initial.csh no_MS '.' Posterior copy 1 ${obs_seq}H

set ms_file   = $1
set local_dir = $2
set kind      = $3
set dim       = $4
set element1  = $5
set yrmoday   = $6

# Import the H* directories information
ls -1 -d H[0-2]* >! hdir_names
set num_times  = `wc -l hdir_names`
set hours      = `cat hdir_names`
rm hdir_names

if (! -d $local_dir) mkdir $local_dir
cd $local_dir

echo "- - - - - - - - - -"
echo "analyses2initial.csh"

if (! -e ${kind}.nc) then
   if (! -e ${kind}_Diag.nc) then
      if (! -e diagnostics.tar) then
         if (! -e diagnostics.tar.gz) then
            echo "starting msrcp at `date`"
            hsi get $ms_file : $ms_file:t 
            # msrcp mss:$ms_file . 
         endif
         echo "starting gunzip at `date`"
         gunzip diagnostics.tar.gz
      endif
      echo "starting tar -x at `date`"
      tar -x -f diagnostics.tar ${kind}_Diag.nc
      if ($status == 0) rm diagnostics.tar &
   endif

   set time = 1
   while ($time <= $num_times[1])
   # Check whether the last file of the 3 has already been created.  If not, generate all 3.
   set out_name = init_analysis_${yrmoday}$hours[$time].nc
   if (! -e $hours[$time]/ice_${out_name}) then
      #-----------
      # CAM
      # Template into which we'll stuff the mean fields: and ensemble average CAM initial file
      if (-e $hours[$time]/caminput_1.nc &&       \
          -e $hours[$time]/clminput_1.nc) then
         cd $hours[$time]
         # clm_ens_avg.f90 needs input.nml for initialization and model version.
         if (! -e input.nml) ln -s ../input.nml .
         pwd
         echo "cam_$out_name clm_$out_name [ice_$out_name]"
         ncra -O -o cam_$out_name caminput*
      else
         echo "analyses2initial.csh; No $hours[$time]/caminput* available; exiting" 
         exit 10
      endif

      # ncks makes the time/record dimension of the out_file to be the largest 
      # of the sizes of the in_ and out_file, if out_file pre-exists.  
      # Prevent this by extracting one timeslot to a new file.
      # (This had -x -v copy in my first successful test, but that's redundant(?) 
      # with the ncwa, where'd I'd like to handle all the 'copy' purging.
      ncks -F -A -d ${dim},${element1} -d time,$time ../${kind}_Diag.nc time_copy_slab.nc 
      if ($status != 0) then
         echo "timeslot $time not found on ../${kind}_Diag.nc"
         exit 20
      endif
      # Exclude the copy dimension and variable from new file, 
      # along with mis-dimensioned P0
      # ncwa will change the rank; use it to remove dimension 'copy'
      ncwa -x -v P0,CopyMetaData -a ${dim} time_copy_slab.nc avgd_copy_out.nc
      # Permute the remaining dimensions into CAM; (time, lev, [s]lat, [s]lon)
      #                                  from DART (time, lat, [s]lon, [s]lev)
      # Requires 3 calls to handle the 3 (staggered) grids
      ncpdq -a lev,lat,lon  avgd_copy_out.nc re-orderedT.nc
      # Is it an FV or eulerian CAM initial file?
      ncdump -v US re-orderedT.nc                               >& /dev/null
      if ($status == 0) then
         # Requires 2 extra calls to handle the staggered grids
         ncpdq -a lev,slat,lon re-orderedT.nc   re-orderedUS.nc
         ncpdq -a lev,lat,slon re-orderedUS.nc  re-ordered.nc
      else
         mv re-orderedT.nc re-ordered.nc
      endif

      # Finally, stuff the contents into the CAM initial file template
      # -A causes variables to be overwritten if there is a name "conflict",
      # without querying the user/script.
      # -x -v $dim removes the *variable* called copy.  (On blueice) it didn't work
      # to try to remove it in the ncwa command.
      ncks -A -x -v ${dim} re-ordered.nc cam_$out_name

      rm time_copy_slab.nc avgd_copy_out.nc re-order* 

      #-----------

      # Added for comparison of forecasts with member 1 vs the ens avg.
      cp clminput_1.nc clm_init_memb1_${yrmoday}$hours[$time].nc 

      # Remove files which won't be removed after archiving in auto_diag2ms_LSF.csh
      rm snow_water_ens.nc input.nml 

      #-----------
      # ICE; ensemble average of the iceinput_#.nc files
      ncea -O -o ice_${out_name} iceinput_[1-9]*.nc
      # also save the first ensemble member for comparison 
      if (-e iceinput_1.nc) cp iceinput_1.nc     ice_init_memb1_${yrmoday}$hours[$time].nc

      cd ..
   endif
   @ time++
   end

endif

echo "- - - - - - - - - -"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


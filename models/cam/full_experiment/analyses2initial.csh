#!/bin/csh

# Extract the mean copy for each state vector field from each of the timeslots 
# on a file (Posterior usually) and put each timeslot into CAM initial file format.
# The non-state fields are ensemble averages of all the caminput_#.nc from the same 
# timeslot.

# The CLM ensemble clminput_#.nc (all non-state fields) is ensemble averaged, 
# and then the snow fields are overwritten using the algorithm from 
# cam#.#.#/models/lnd/clm2/src/main/snowdp2lev.F90
# packaged as an NCO script, which is used by ncap2.
# ncap2 may need to be updated to NCO 3.9.4 or later.

# set echo verbose


# Called from auto_diag2ms_LSF.csh with
#   ../../analyses2initial.csh no_MS '.' Posterior copy 1 ${obs_seq}H              >>& $saved

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
   while ($time <= $num_times[1])
      set out_name = init_analysis_${yrmoday}$hours[$time].nc
      #-----------
      # CAM
      # Template into which we'll stuff the mean fields: and ensemble average CAM initial file
      if (-e $hours[$time]/caminput_1.nc &&       \
          -e $hours[$time]/clminput_1.nc) then
         cd $hours[$time]
         # clm_ens_avg.f90 needs input.nml for initialization and model version.
         if (! -e input.nml) ln -s ../input.nml .
         pwd
         echo cam_$out_name clm_$out_name
         ncra -O -o cam_$out_name caminput*
      else
         echo "analyses2initial.csh; No $hours[$time]/caminput* available; exiting" 
         exit
      endif

      # ncks makes the time/record dimension of the out_file to be the largest 
      # of the sizes of the in_ and out_file, if out_file pre-exists.  
      # Prevent this by extracting one timeslot to a new file.
      # (This had -x -v copy in my first successful test, but that's redundant(?) 
      # with the ncwa, where'd I'd like to handle all the 'copy' purging.
      ncks -F -A -d ${dim},${element1} -d time,$time ../${kind}_Diag.nc time_copy_slab.nc 
      if ($status != 0) then
         echo "timeslot $time not found on ../${kind}_Diag.nc"
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
      # CLM; ensemble average of the clminput files.   ncra can't be used because files
      #      have no record dimension.  This will have incorrect snow and water fields.
      ncea -O -o clm_ens_avg.nc clminput_[1-9]*.nc 

      # Create a file of the ensemble of the snow and water fields which need to be fixed.
      set cat_flds = 'DZSNO,H2OSNO,H2OSOI_LIQ,H2OSOI_ICE,SNLSNO,SNOWDP,T_SOISNO'
      set cat_flds = "${cat_flds},snw_rds,qflx_snofrz_lyr,mss_bcpho,mss_bcphi,mss_ocpho,mss_ocphi"
      set cat_flds = "${cat_flds},mss_dst1,mss_dst2,mss_dst3,mss_dst4"
      set cat_flds = "${cat_flds},flx_absdv,flx_absdn,flx_absiv,flx_absin"
      ncecat -u ensemble -v ${cat_flds} -o snow_water_ens.nc clminput_[1-9]*.nc

      # Fix the snow fields with fortran program clm_ens_avg.f90 (different versions for 
      # different CLM initial file versions.)
      # This reads snow_water_ens.nc and writes into clm_ens_avg.nc.
      ../../../clm_ens_avg 

      # ncea prunes unused dimensions from the output file, but CLM needs them,
      # so insert all the good, averaged, snow and water fields into a copy of an ensemble member, 
      # which has a complete header.  Also exclude the timemgr fields.
      cp clminput_1.nc clm_${out_name}
      ncks -A -x -v '^timemgr' -o clm_${out_name} clm_ens_avg.nc

# Remove files which won't be removed after archiving in auto_diag2ms_LSF.csh
      rm snow_water_ens.nc input.nml dart*

      #-----------
      # ICE; no averaging available yet; use the first ensemble member as the "analysis"
      cp iceinput_1.tar     ice_init_memb1_${yrmoday}$hours[$time].tar

      cd ..
      @ time++
   end

endif

echo "- - - - - - - - - -"
exit

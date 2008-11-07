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

# Upgrade; can these variables be filled by looking at what H* directories exist?
# Then they wouldn't have to be hardwired here.
set num_times  = 4
set hours      = (06 12 18 24)

if (! -d $local_dir) mkdir $local_dir
cd $local_dir

if (! -e ../../snow.nco) then
   echo '/* '                                                                    >! snow.nco
   echo 'execute via'                                                            >> snow.nco 
   echo 'ncap2 -O -S snow.nco -o clminput_avg_goodsno.nc clminput_avg_badsno.nc' >> snow.nco
   echo ' '                                                                      >> snow.nco
   echo 'where clminput_avg_badsno.nc is the ensemble average '                  >> snow.nco

   echo 'of the clminput_#.nc files'                                             >> snow.nco
   echo '---------------------------'                                            >> snow.nco
   echo '*/ '                                                                    >> snow.nco
   echo ' '                                                                      >> snow.nco
   echo '*cols = $column.size  ;'                                                >> snow.nco
   echo '*landkind = cols1d_ityplun ;'                                           >> snow.nco
   echo '*SNOWDPr = SNOWDP ;'                                                    >> snow.nco
   echo ' '                                                                      >> snow.nco
   echo '*SNLSNOr[column]        = 0    ;'                                       >> snow.nco
   echo '*DZSNOr [column,levsno] = 0.0d ;'                                       >> snow.nco
   echo '*ZSNOr  [column,levsno] = 0.0d ;'                                       >> snow.nco
   echo '*ZISNOr [column,levsno] = 0.0d ;'                                       >> snow.nco
   echo ' '                                                                      >> snow.nco
   echo 'for (*c=0s ; c< cols; c++) {'                                           >> snow.nco
   echo '   if ((SNOWDPr(c)>= 0.01d) && (SNOWDPr(c)<= 0.03d)) {'                 >> snow.nco
   echo '      SNLSNOr(c) = -1 ;'                                                >> snow.nco
   echo '      DZSNOr(c, 4)  = SNOWDPr(c);'                                      >> snow.nco
   echo '    } else if ((SNOWDPr(c)> 0.03d) && (SNOWDPr(c)<= 0.04d)) {'          >> snow.nco
   echo '      SNLSNOr(c) = -2 ;'                                                >> snow.nco
   echo '      DZSNOr(c, 3) = SNOWDPr(c)/2.0d ;'                                 >> snow.nco
   echo '      DZSNOr(c, 4) = DZSNOr(c, 3) ;'                                    >> snow.nco
   echo '    } else if ((SNOWDPr(c)> 0.04d) && (SNOWDPr(c)<= 0.07d)) {'          >> snow.nco
   echo '      SNLSNOr(c) = -2 ;'                                                >> snow.nco
   echo '      DZSNOr(c, 3) = 0.02d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 4) = SNOWDPr(c)- DZSNOr(c, 3) ;'                        >> snow.nco
   echo '    } else if ((SNOWDPr(c)> 0.07d) && (SNOWDPr(c)<= 0.12d)) {'          >> snow.nco
   echo '      SNLSNOr(c) = -3 ;'                                                >> snow.nco
   echo '      DZSNOr(c, 2) = 0.02d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 3) = (SNOWDPr(c)- DZSNOr(c, 2))/2.0d ;'                 >> snow.nco
   echo '      DZSNOr(c, 4) = DZSNOr(c, 3) ;'                                    >> snow.nco
   echo '    } else if ((SNOWDPr(c)> 0.12d) && (SNOWDPr(c)<= 0.18d)) {'          >> snow.nco
   echo '      SNLSNOr(c) = -3 ;'                                                >> snow.nco
   echo '      DZSNOr(c, 2) = 0.02d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 3) = 0.05d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 4) = SNOWDPr(c)- DZSNOr(c, 2) - DZSNOr(c, 3) ;'         >> snow.nco
   echo '    } else if ((SNOWDPr(c)> 0.18d) && (SNOWDPr(c)<= 0.29d)) {'          >> snow.nco
   echo '      SNLSNOr(c) = -4 ;'                                                >> snow.nco
   echo '      DZSNOr(c, 1) = 0.02d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 2) = 0.05d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 3) = (SNOWDPr(c)- DZSNOr(c, 1) - DZSNOr(c, 2))/2.0d ;'  >> snow.nco
   echo '      DZSNOr(c, 4) = DZSNOr(c, 3) ;'                                    >> snow.nco
   echo '    } else if ((SNOWDPr(c)> 0.29d) && (SNOWDPr(c)<= 0.41d)) {'          >> snow.nco
   echo '      SNLSNOr(c) = -4 ;'                                                >> snow.nco
   echo '      DZSNOr(c, 1) = 0.02d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 2) = 0.05d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 3) = 0.11d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 4) = SNOWDPr(c)- DZSNOr(c, 1) - DZSNOr(c, 2) - DZSNOr(c, 3) ;' >> snow.nco
   echo '    } else if ((SNOWDPr(c)> 0.41d) && (SNOWDPr(c)<= 0.64d)) {'          >> snow.nco
   echo '      SNLSNOr(c) = -5 ;'                                              >> snow.nco
   echo '      DZSNOr(c, 0) = 0.02d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 1) = 0.05d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 2) = 0.11d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 3) = (SNOWDPr(c)- DZSNOr(c, 0) - DZSNOr(c, 1) - DZSNOr(c, 2))/2.0d ;' >> snow.nco
   echo '      DZSNOr(c, 4) = DZSNOr(c, 3) ;'                                    >> snow.nco
   echo '    } else if (SNOWDPr(c)> 0.64d) {'                                    >> snow.nco
   echo '      SNLSNOr(c) = -5 ;'                                                >> snow.nco
   echo '      DZSNOr(c, 0) = 0.02d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 1) = 0.05d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 2) = 0.11d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 3) = 0.23d ;'                                           >> snow.nco
   echo '      DZSNOr(c, 4) = SNOWDPr(c)- DZSNOr(c, 0)- DZSNOr(c, 1)- DZSNOr(c, 2)- DZSNOr(c, 3) ;' >> snow.nco
   echo '    }'                                                                  >> snow.nco
   echo '    '                                                                   >> snow.nco
   echo '//   lake points do not get any snow distributed in snowdp2lev.'        >> snow.nco 
   echo '   if ( landkind(c) != 3) {'                                            >> snow.nco
   echo '      for ( *j = 4; j>5+SNLSNOr(c); j--) {'                             >> snow.nco
   echo '         ZSNOr(c,j)    = ZISNOr(c,j) - DZSNOr(c,j) * 0.5d ;'            >> snow.nco
   echo '         ZISNOr(c,j-1) = ZISNOr(c,j) - DZSNOr(c,j) ;'                   >> snow.nco
   echo '      }'                                                                >> snow.nco
   echo '      // j has been decremented to 0'                                   >> snow.nco
   echo '      // This was taken out of the (shortened) loop to avoid referencing ZISNOr(c,-1)' >> snow.nco
   echo '      ZSNOr(c,j)    = ZISNOr(c,j) - DZSNOr(c,j) * 0.5d ;'               >> snow.nco
   echo '   }'                                                                   >> snow.nco
   echo '}'                                                                      >> snow.nco
   echo ' '                                                                      >> snow.nco
   echo 'SNLSNO = SNLSNOr;'                                                      >> snow.nco
   echo 'DZSNO  = DZSNOr ;'                                                      >> snow.nco
   echo 'ZSNO   = ZSNOr  ;'                                                      >> snow.nco
   echo 'ZISNO  = ZISNOr ;'                                                      >> snow.nco

   mv snow.nco ../..
endif

echo "- - - - - - - - - -"
echo "analyses2initial.csh":

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
      set out_name = init_analysis_${yrmoday}$hours[$time].nc
      # Template into which we'll stuff the mean fields
#     Want ensemble average CAM initial file
      if (-e H$hours[$time]/caminput_1.nc &&       \
          -e H$hours[$time]/clminput_1.nc) then
         cd H$hours[$time]
         pwd
         echo cam_$out_name clm_$out_name
         ncra -o cam_$out_name caminput*
      else
         echo "analyses2initial.csh; No H$hours[$time]/caminput* available; exiting" 
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
      ncpdq -a lev,slat,lon re-orderedT.nc   re-orderedUS.nc
      ncpdq -a lev,lat,slon re-orderedUS.nc  re-orderedVS.nc
      # Finally, stuff the contents into the CAM initial file template
      # -A causes variables to be overwritten if there is a name "conflict",
      # without querying the user/script.
      # -x -v $dim removes the *variable* called copy.  (On blueice) it didn't work
      # to try to remove it in the ncwa command.
      ncks -A -x -v ${dim} re-orderedVS.nc cam_$out_name

      # CLM; ensemble average of the clminput files, 
      ncea -O -o clminput_avg_badsno.nc clminput_*.nc 
      # Fix the snow fields with NCO script based on algorithm in
      # Cam3/cam3.5/models/lnd/clm2/src/main/snowdp2lev.F90
      ncap2 -O -S ../../../snow.nco -o clm_$out_name clminput_avg_badsno.nc

      rm time_copy_slab.nc avgd_copy_out.nc re-order* clminput_avg_badsno.nc
      cd ..
      @ time++
   end

endif

echo "- - - - - - - - - -"
exit

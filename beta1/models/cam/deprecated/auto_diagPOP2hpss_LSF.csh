#!/bin/tcsh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# script for copying diagnostics files to mass store.

#### LSF options for BSUB
### -J      job name    (master script job.csh presumes filter.xxxx.log)
### -o      output listing filename
### -P      account number
### -q # Queue name    regular   economy  standby     long   
#        proclim            32        16        8        2      
#        timelim         6 hrs        18       48   5 days  
#        # jobs/person       2         -        -        2
### -n      number of tasks (processors)
### -x      exclusive use of node
### -R "span[ptile=(num procs you want on each node)]"
#
#BSUB -J auto_diag2hpss
#BSUB -o auto_diag2hpss.%J.log
#BSUB -e auto_diag2hpss.%J.err
#BSUB -P xxxxxxxx
#BSUB -q share
#BSUB -W 2:00
#BSUB -n 1
#BSUB -R "span[ptile=1]"


set compress = true

set diag_name = diagnostics.tar
set saved = saved_diagnostics

if ($?LS_SUBCWD) then
   cd $LS_SUBCWD
endif

touch $saved
echo '------------------------------------------------------------' >> $saved
echo 'auto_diag2hpss_LSF starts in' >> $saved
pwd                               >> $saved
date                              >> $saved

set direct = `pwd`
set obs_seq = $direct:t

cd ..
set direct = `pwd`
set case = $direct:t

cd ..
set direct = `pwd`
set exp_dir = $direct:t

cd $case/${obs_seq}
set hpss_dir = /RAEDER/DAI/${exp_dir}/$case/${obs_seq}

# make sure this exists before we start
hsi mkdir ${hpss_dir}

# IBM tar requires 1 entry/line in list of things to exclude from tar
echo DART                >! tar_excl_list
echo CAM                 >> tar_excl_list
echo CLM                 >> tar_excl_list
echo ICE                 >> tar_excl_list
# batch* are the files into which DART,CAM,CLM are archived by auto_re2hpss_LSF.csh,
# which is run at the same time as this script.
echo 'batch*'            >> tar_excl_list
echo $saved              >> tar_excl_list

## Added to make mean easily accessible in the form of a CAM initial file
echo 'cam_analyses.tar'  >> tar_excl_list
echo 'H[0-9]*'           >> tar_excl_list
echo 'H_*'               >> tar_excl_list

#-----------------------------
# Stuff the Posterior mean fields into CAM initial files.
# Tar together with a CLM initial file.
# Arguments to analyses2initial.csh are
#   set hpss_file   = $1   script searches for local Posterior_Diag.nc first, so give a dummy MS name.
#   set local_dir = $2
#   set kind      = $3
#   set dim       = $4
#   set element1  = $5
#   set yrmoday   = $6   set to $obs_seq instead of yyyymmdd, since obs_seq is easily available

# Save out history files from H* directories before calculating CAM/CLM analyses and deleting H*
# compressing saves a factor of 12.
# Can't do with first CAM 3.6.71; empty_htapes conflicted with ENDOFRUN
 ls H[012]*/*.h0.* >& /dev/null
if ($status == 0) then
   gzip H[012]*/*.h0.*
   tar -c -f H_all.h0.gz.tar H[012]*/*.h0.*
   set ar_status = $status
   if ($ar_status == 0 && -e H_all.h0.gz.tar) then
      hsi put H_all.h0.gz.tar : ${hpss_dir}/H_all.h0.gz.tar                    >>& $saved  &
      set ar_status = $status
      if ($ar_status == 0) rm H[012]*/*.h0.* 
   endif
   if ($ar_status != 0) echo 'ARCHIVING of H[012]*.h0.gz.tar FAILED' >> & $saved
else
    echo 'ARCHIVING of H[012]*.h0.gz.tar FAILED; no files available' >> & $saved
endif 

# Save out coupler history files from H* directories before calculating CAM/CLM analyses and deleting H*
# Compressing saves a factor of ??.
# File names have the form H##/FV_2deg_greg-O2-POP1-1.cpl.ha2x6h[r].2006-12-03.nc   ## = 06,12,18,24
# These fields were identified from the sample input file provided by Yeager and/or Lindsey. 


# ? Keep time_bnds too

set flds = 'time,doma_lat,doma_lon,doma_area,doma_mask'
set flds = "${flds},a2x6h_Faxa_swndr,a2x6h_Faxa_swvdr,a2x6h_Faxa_swndf,a2x6h_Faxa_swvdf"
set flds = "${flds},a2x6h_Faxa_rainc,a2x6h_Faxa_rainl,a2x6h_Faxa_snowc,a2x6h_Faxa_snowl"
set flds = "${flds},a2x6h_Sa_z,a2x6h_Sa_u,a2x6h_Sa_v,a2x6h_Sa_tbot,a2x6h_Sa_ptem,a2x6h_Sa_shum"
set flds = "${flds},a2x6h_Sa_pbot,a2x6h_Sa_dens,a2x6h_Sa_pslv,a2x6h_Faxa_lwdn"

set memb = 1
set more = true
while ($more == true )
   ls H[012]*/*-${memb}.cpl.ha2x6h*.* >& /dev/null
   if ($status != 0) then
      set more = false
      if ($memb == 1) then
         echo 'ARCHIVING of H[012]/*cpl.ha2x6h.* FAILED; no files available' >> & $saved
      endif
   else
      # Extract filenames' pieces
      set hists = `ls H[012]*/*-${memb}\.cpl\.ha2x6h*`
      set date = $hists[1]:r:e
      set case = $hists[1]:t:r:r:r:r

      # Ensemble average of the 4 6-hour averages for each ensemble member.


      ncra -O -v ${flds} -o ${case}.cpl.ha2x1davg.${date}.nc $hists

      echo 'ncra created average file:'                                     >> & $saved
      ls ${case}.cpl.ha2x1davg.${date}.nc                                   >> & $saved
      echo 'from'                                                           >> & $saved
      ls -l $hists                                                          >> & $saved

      # ncea prunes unused dimensions (here; time) from the output file.  
      # That's OK because lots of these averaged
      # files will be concatenated into the single file which datm/cplr wants to see.
      # The time dimension can be reinstated then.
      # So don't remove the averaged files after they've been archived; 
      # concatenate them together into months while bundling the obs_seq.final files.
      ncap2 -O -s 'time[$time]={.5}' ${case}.cpl.ha2x1davg.${date}.nc ${case}.cpl.ha2x1davg.${date}.nc 

      # Concatenate the times for each day into a single file using record (time) concatenator.
      ncrcat -v ${flds} -o ${case}.cpl.ha2x1dx6h.${date}.nc $hists
      echo 'ncrcat created daily time series file:'                                     >> & $saved
      ls ${case}.cpl.ha2x1dx6h.${date}.nc                                   >> & $saved
      echo 'from'                                                           >> & $saved
      ls -l $hists                                                          >> & $saved

      # The times on the concatenated cplr history file are time = 0.25, 0.25, 0.25, 0.25 ;
      #   as measured from    time:units = "days since 2006-12-02 00:00:00" ;
      # ? Fix that when generating, or 
      #   Fix for each day here, then add day information during final concatenation into monthly files.
      ncap2 -O -s 'time[$time]={0.25,.5,.75,1.0}' ${case}.cpl.ha2x1dx6h.${date}.nc \
            ${case}.cpl.ha2x1dx6h.${date}.nc

      # Compress the daily avg and time series files for this ens member while working on the next
      gzip ${case}.cpl.ha2x1d* &

   endif
   @ memb++
end

wait
tar -c -f H_cplr.ha2x1dx6h.gz.tar *.ha2x1dx6h*gz &
tar -c -f H_cplr.ha2x1davg.gz.tar *.ha2x1davg*gz &
wait
tar -c -f H_cplr.ha2x1d.gz.tar *.ha2x1d[ax]*tar
set ar_status = $status
if ($ar_status == 0 && -e H_cplr.ha2x1d.gz.tar) then
   hsi put H_cplr.ha2x1d.gz.tar : ${hpss_dir}/H_cplr.ha2x1d.gz.tar                    >>& $saved  
   set ar_status = $status
   # Leave ha2x2davg.gz.tar here for monthly archiving, like obs_seq.final
   if ($ar_status == 0) then
      rm H[012]*/*.ha2x6h*  *.ha2x1dx6h*.gz *.ha2x1davg*.gz &
      echo "SUCCEEDED archiving H_cplr.ha2x1d.gz.tar" >>& $saved
   endif
endif
if ($ar_status != 0) echo 'ARCHIVING of *.ha2x1d.gz.tar FAILED' >> & $saved

#-------------------------------------
# Archive the analyses.
if (-e ../../analyses2initial.csh) then
   # analyses2initial.csh needs CAM initial files to average and receive the analyses 
   # from Posterior_Diag.nc.
   # They should have been saved during the assimilation and be living in the 
   # exp/obs_####/H## directories.
   ls H*/caminput_1.nc > /dev/null
   set stat = $status
   if ($stat == 0) then
      ls H*/clminput_1.nc > /dev/null
      set stat = $status
   endif
   if ($stat != 0) then
      echo "H*/c[al]minput_* not available"                           >>& $saved
      echo "H*/c[al]minput_* not available" >&  ANALYSES_NOT_SAVED
   else
#   if ($status == 0) ls H*/clminput_1.nc > /dev/null
#   if ($status != 0) then
#      echo "H*/c[al]minput_* not available; exiting"                           >>& $saved
#      echo "H*/c[al]minput_* not available; exiting" >&  ANALYSES_NOT_SAVED
#      exit
#   endif
   
      set num_anal = `ls H[0-2]*/cam_init_*`
      set tar_stat = 0
      if (! -e cam_analyses.tar) then
         if ($#num_anal < 4) then
            echo " "                                                                >>& $saved
            ../../analyses2initial.csh no_MS '.' Posterior copy 1 ${obs_seq}        >>& $saved
            set num_anal = `ls H[0-2]*/cam_init_*`
         endif
   
         tar -c -f cam_analyses.tar H[0-2]*/{c,ice}*_init_* 
         set tar_stat = $status
      endif
      # This section requires that old/failed save_diagnostic files be left in place.
      # cam_init_analysis is part of a filename constructed in analyses2initial.csh
      # and printed for each H# directory.
# Sometimes H06 doesn't get written to saved_diagnostics by analyses2initial.csh
#      grep cam_init_analysis saved_diagnostics >&! cam_init
#      set num_anal = `wc -l cam_init`
# This assumes 4 analysis times/obs_seq file
#      if ($tar_stat == 0 && $num_anal[1] > 3)  \
      if ($tar_stat == 0 )  \
         hsi put cam_analyses.tar : ${hpss_dir}/cam_analyses.tar                    >>& $saved  
      set list = `ls -l cam_analyses.tar`
      set local_size = $list[5]
      set list = `msls -l ${hpss_dir}/cam_analyses.tar`
      set hpss_size = $list[5]
      echo " cam_analyses.tar local_size = $local_size, hpss_size = $hpss_size"       >> $saved
   
      if ($local_size == $hpss_size) then
         echo "Archived $hpss_dir/cam_analyses.tar                                " >> $saved
         echo '    REMOVING H[0-9]* and cam_analyses.tar '                        >> $saved
         rm -rf H[0-9]*/[ci]* cam_analyses.tar cam_init
      else
         echo "hsi put of ${hpss_dir}/cam_analyses.tar  failed; "                     >> $saved
         echo 'NOT removing H[0-9]* and cam_analyses.tar '                        >> $saved
      endif
   endif
else
   echo "NO analyses2initial.csh, so no CAM initial file format analyses created"
endif

#-----------------------------

if (! -e $diag_name && ! -e ${diag_name}.gz) tar -c -v -f $diag_name -X tar_excl_list * >>& $saved
rm tar_excl_list

if ($compress == true) then
   if (-e $diag_name) then
      gzip $diag_name                         >>& $saved
      set diag_name = ${diag_name}.gz
   else if (-e ${diag_name}.gz) then
      set diag_name = ${diag_name}.gz
   else
      echo "$diag_name does not exist at gzip" >> $saved
      exit
   endif
endif

echo "files will be written to ${hpss_dir}/${diag_name}" >> $saved

hsi put ${diag_name} : ${hpss_dir}/${diag_name} >>& $saved

set list = `ls -l $diag_name`
set local_size = $list[5]
set list = `msls -l ${hpss_dir}/${diag_name}`
set hpss_size = $list[5]
echo " ${diag_name} local_size = $local_size, hpss_size = $hpss_size" >> $saved

if ($local_size == $hpss_size) then
   echo "hsi put of $hpss_dir/$diag_name succeeded; REMOVING $diag_name and P*Diag.nc " >> $saved
   rm $diag_name P*Diag.nc
else
   echo "hsi put of ${hpss_dir}/$obs_seq  failed; " >> $saved
   echo "NOT removing $diag_name and P*.nc"      >> $saved
endif

chmod 444 $saved

wait
if ($ar_status == 0) rm H_all.h0.gz.tar  *.ha2x1dx6h*.gz.tar

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


#!/bin/bash
###############################################################
# Author: Jingjing Liang, on loan to UT Austin 
#
# This runs NOAH-MP 3.6 and the filter on the UT 'Lonestar5' machine.
#
##############################################################

#  change here to initialize the Noah-MP DA

startdate='2003-08-12 00:00:00'
enddate='2004-01-30 00:00:00'
interval=86400                 # seconds,which means a day
numens=40                      # ensemble size -- used during perturbation

obsdir='/work/04270/l_jing/lonestar/DART_CAM/merge/output'
subout='outfiles'

noahdir='/scratch/04270/l_jing/DART_NoahMP'
outdir="/scratch/04270/l_jing/DART_NoahMP/output/${subout}"
DARTDIR="/home1/04270/l_jing/DART_NoahMP"
ana_pre_dir='/scratch/04270/l_jing/DART_NoahMP/output/outfiles/ana_pre_MOD_GRA'
inf_dir='/scratch/04270/l_jing/DART_NoahMP/output/outfiles/inf_MOD_GRA'
subhist='/hist_MOD_GRA'
subrest='/rest_MOD_GRA'
sublogs='/logs_MOD_GRA'

################################################################ 

echo "Begin Noah-MP DA    ctime = `date`"
echo ""

ddays=`expr '(' $(date +%s -d "$enddate") - $(date +%s -d "$startdate") ')' / 86400 + 1`
dtemp=0
while (($dtemp<$ddays))
do
   #----------------------------------------------------
   # retrieve model time info 
   cds=`expr $dtemp \* $interval`
   cyy=`date -d "$startdate $cds seconds" +"%Y"`
   cmm=`date -d "$startdate $cds seconds" +"%m"`
   cdd=`date -d "$startdate $cds seconds" +"%d"`
  # cHH=`date -d "$startdate $cds seconds" +"%H"`
   cHH='00'
   echo " " 
   echo "Start DA cycle on ${cyy}-${cmm}-${cdd}:${cHH}" 

   #----------------------------------------------------
   #  update date
   dtemp=`expr $dtemp + 1`
   nds=`expr $dtemp \* $interval`
   nyy=`date -d "$startdate $nds seconds" +"%Y"`
   nmm=`date -d "$startdate $nds seconds" +"%m"`
   ndd=`date -d "$startdate $nds seconds" +"%d"`
  # nHH=`date -d "$startdate $nds seconds" +"%H"` 
   nHH='00'
   echo "Start DA cycle on ${nyy}-${nmm}-${ndd}:${nHH}" 
   #----------------------------------------------------
   # Noah-MP namelist modify and implementation
   for(( inst = 1; inst <= ${numens}; inst++ ))
   do
   {  
      inst_str=`printf %02d ${inst}`
     # if [[ $cds -lt 86400 ]]; then
      if [[ $cds -lt 0 ]]; then
        cd ${noahdir}/hrldas_${inst_str}
        inifile1="RESTART.${cyy}${cmm}${cdd}${cHH}_DOMAIN1"
sed -i "14s:.*:RESTART_FILENAME_REQUESTED = '/scratch/04270/l_jing/DART_NoahMP/hrldas_${inst_str}/output/${inifile1}':" namelist.hrldas
      else
        cd ${noahdir}/hrldas_${inst_str}
        inifile="RESTART.${cyy}${cmm}${cdd}${cHH}_DOMAIN1_${inst_str}"
sed -i "1,\$s#RESTART_FILENAME_REQUESTED.*\$#RESTART_FILENAME_REQUESTED  = \"${outdir}\/${subrest}\/${inifile}\"#g" namelist.hrldas
      fi
      sed -i "1,\$s#START_YEAR.*\$#START_YEAR = $cyy #g" namelist.hrldas
      sed -i "1,\$s#START_MONTH.*\$#START_MONTH = $cmm #g" namelist.hrldas
      sed -i "1,\$s#START_DAY.*\$#START_DAY = $cdd #g" namelist.hrldas
      sed -i "1,\$s#START_HOUR.*\$#START_HOUR = $cHH #g" namelist.hrldas
      ./noahmp_hrldas.exe >noahmp.log_$cyy$cmm$cdd 
   }&
   done

#-------------------------------------------------------------------------------
# this part could do the task through array, but needs more time for submit&wait
# if submit this main.sh as a job, the sentance with LJJ should be used 
   #   cd /home1/04270/l_jing/DART_NoahMP/launcher
   #   echo "*****RUNNING******* " > lock
   #   sbatch launcher.slurm
   #   count=1
   #   for (( count=1;count<=100;count++ )) do
   #       if [ -f lock ]
   #       then
   #           echo "    WAIT PREVIOUS JOB TO FINISH... "
   #           sleep 10
   #       else
   #           break
   #       fi
   #    done
   #
   #echo ''
   echo "     move output files ..."
   sleep 15
#-------------------------------------------------------------------------------
   cd $noahdir
   inst=1
   while ((inst<=numens))
   do
      inst_str=`printf %02d ${inst}`
      cd ${noahdir}/hrldas_${inst_str}/output_tws
      restfile="RESTART.${nyy}${nmm}${ndd}${nHH}_DOMAIN1"
      histfile="${cyy}${cmm}${cdd}${cHH}.LDASOUT_DOMAIN1"
      mv ${histfile} "${outdir}${subhist}/${histfile}_${inst_str}"
      mv ${restfile} "${outdir}${subrest}/${restfile}_${inst_str}"
      cd ..
      mv "noahmp.log_$cyy$cmm$cdd" "${outdir}${sublogs}/noahmp.log_${cyy}${cmm}${cdd}_${inst_str}"
      inst=`expr $inst + 1`
   done
 
   echo "  Finish Noah-MP on ${nyy}-${nmm}-${ndd}:${nHH}"     
   #----------------------------------------------------
   # check if observations are available for this day
    obsfilename="$obsdir/MOD_GRA_obs_seq.${nyy}${nmm}${ndd}" #it's for MOD_GRA 
   if [ ! -e $obsfilename ]
   then 
      echo "    No obs on ${nyy}-${nmm}-${ndd} - continue without DA"
      continue
   fi
   #---------------------------------------------------
   # to do the filter
   # First, inflation
   # IF we are doing inflation, we must take the output inflation files from
   # the previous cycle and rename them for input to the current cycle.
   
   # Should the setup script just create input inflation files so we don't
   # have to screw with changing the namelist after the first execution
   # (which traditionally reads from the namelist, not the file)
   
   # If the file exists, just link to the new expected name.
   # the expected names have a _d0? inserted before the file extension
   # if there are multiple domains.
   # If the file does not exist, filter will die and issue a very explicit
   # death message.
   cd $DARTDIR
   rm -f input_priorinf_mean.nc input_priorinf_sd.nc
 
    # Checking for a prior inflation mean file from the previous assimilation.
 
    (ls -rt1 ${inf_dir}/analysis_priorinf_mean.*.nc | tail -n 1 > latestfile) 
    nfiles=`cat latestfile | wc -l`
 
    if [[ $nfiles -gt 0 ]]; then
       latest=`cat latestfile`
       ln -vs $latest input_priorinf_mean.nc
    fi
 
    # Checking for a prior inflation sd file from the previous assimilation.
 
    (ls -rt1 ${inf_dir}/analysis_priorinf_sd.*.nc | tail -n 1 > latestfile) 
    nfiles=`cat latestfile | wc -l`
 
    if [[ $nfiles -gt 0 ]]; then
       latest=`cat latestfile`
       ln -vs $latest input_priorinf_sd.nc
    fi
 
   # Then, begin to do the filter
   # write the name of restart file to the input_file_list/output_file_list
   ls ${outdir}${subrest}/RESTART.${nyy}${nmm}${ndd}${nHH}_DOMAIN1* > $DARTDIR/input_file_list.txt         
   cd $DARTDIR
   cp input_file_list.txt output_file_list.txt

   rm -f obs_seq.out
   ln -s ${obsdir}/MOD_GRA_obs_seq.$nyy$nmm$ndd obs_seq.out
echo 'do filter'
   ./filter || exit 1
echo 'finish filter'
   mv preassim_priorinf_mean.nc ${inf_dir}/preassim_priorinf_mean.${nyy}${nmm}${ndd}${nHH}.nc
   mv preassim_priorinf_sd.nc   ${inf_dir}/preassim_priorinf_sd.${nyy}${nmm}${ndd}${nHH}.nc
   mv analysis_priorinf_mean.nc ${inf_dir}/analysis_priorinf_mean.${nyy}${nmm}${ndd}${nHH}.nc
   mv analysis_priorinf_sd.nc   ${inf_dir}/analysis_priorinf_sd.${nyy}${nmm}${ndd}${nHH}.nc
   mv obs_seq.final           /scratch/04270/l_jing/DART_NoahMP/observation/MOD_GRA/obs_seq.final_${nyy}${nmm}${ndd}
   #to run the ./dart_to_noah, use it when assimilate MODIS
   for(( inst = 1; inst <= ${numens}; inst++ ))
   do
   {
      inst_str=`printf %04d ${inst}`
      res_str=`printf %02d ${inst}`
      mv analysis_member_${inst_str}.nc ${ana_pre_dir}/analysis_member_${inst_str}_${nyy}${nmm}${ndd}${nHH}.nc
      mv preassim_member_${inst_str}.nc ${ana_pre_dir}/preassim_member_${inst_str}_${nyy}${nmm}${ndd}${nHH}.nc
      rm -f analysis_input.nc
      ln -s ${ana_pre_dir}/analysis_member_${inst_str}_${nyy}${nmm}${ndd}${nHH}.nc analysis_input.nc
      rm -f restart_input.nc
      ln -s ${outdir}${subrest}/RESTART.${nyy}${nmm}${ndd}${nHH}_DOMAIN1_${res_str} restart_input.nc
      ./dart_to_noah > dart_to.log
}  
   done
done

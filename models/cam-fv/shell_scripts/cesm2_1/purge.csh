#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# $Id:$

#==========================================================================
# This script removes assimilation output from $DOUT_S_ROOT ($local_archive) after:
# 1) it has been repackaged by repack_st_archive.csh,
# 2) globus has finished moving files to data_campaign storage.,
# 3) other output has been left behind by all of those processes (i.e in rundir).
# >>> Check all final archive locations for complete sets of files
#     before running this.  <<<
#     This can be done by using ./pre_purge.csh, 
#     or see "Look in" below.
# Submit from $CASEROOT.

#-----------------------------------------

if (! -f CaseStatus) then
   echo "ERROR: this script must be run from the CESM CASEROOT directory"
   exit 1
endif

#--------------------------------------------
# Default values of which file sets to purge.
set do_forcing     = 'true'
set do_restarts    = 'true'
set do_history     = 'true'
set do_state_space = 'true'
set do_rundir      = 'true'

# "components" = generic pieces of CESM (used in the archive directory names).
# "models" = component instance names (models, used in file names).
set components     = (lnd  atm ice  rof)
set models         = (clm2 cam cice mosart)


if ($#argv != 0) then
   # Request for help; any argument will do.
   echo "Usage:  "
   echo "Before running this script"
   echo "    Run repack_st_archive.csh. "
   echo "    Confirm that it put all the output where it belongs using pre_purge.csh. "
   echo "    Edit the script to be sure that your desired "
   
   echo "    file types and model/components will be purged"
   echo "Call by user or script:"
   echo "   purge.csh data_proj_space_dir data_campaign_dir year mo [do_this=false] ... "
   exit
endif

source ./data_scripts.csh

set yr_mo = `printf %4d-%02d ${data_year} ${data_month}`
set local_arch     = `./xmlquery DOUT_S_ROOT --value`

# These arguments to turn off parts of the archiving must have the form do_obs_space=false ;
# no quotes, no spaces.
#    if ($?5) set $5
#    if ($?6) set $6
#    if ($?7) set $7

#--------------------------------------------

set lists = logs/rm_${yr_mo}.lists
if (-f ${lists}.gz) mv ${lists}.gz ${lists}.gz.$$
pwd > $lists

#------------------------
# Look in ${data_proj_space}/${CASE}/cpl/hist/$INST
# e.g. ${pr}/${casename}/cpl/hist/0080
if ($do_forcing == true) then
   cd ${local_arch}/cpl/hist
   pwd >>& ${local_arch}/$lists

   # The original cpl hist (forcing) files,
   # which have been repackaged into $data_proj_space.
   foreach type (ha2x3h ha2x1h ha2x1hi ha2x1d hr2x)
      echo "Forcing $type" >>& ${local_arch}/$lists
      ls ${data_CASE}.cpl_*.${type}.${yr_mo}-*.nc >>& ${local_arch}/$lists
      rm ${data_CASE}.cpl_*.${type}.${yr_mo}-*.nc >>& ${local_arch}/$lists
   end

   # These files were moved from $data_proj_space storage to be input to ncrcat.
   # (2020-1-18; that was before implementing Nancy's suggestion to leave
   #  the yearly file on $pr.)
   # An updated version was put into $data_proj_space, so the locals are not needed.
#    echo "Forcing $year" >>& ${local_arch}/$lists
#    ls *${year}.nc >>& ${local_arch}/$lists
#    rm *${year}.nc >>& ${local_arch}/$lists

   echo "Forcing \*.eo" >>& ${local_arch}/$lists
   ls *.eo >>& ${local_arch}/$lists
   rm *.eo >>& ${local_arch}/$lists

   cd ${local_arch}
endif

#------------------------
# Look in ${data_campaign}/${data_CASE}/rest/$yr_mo
# E.g. ${data_campaign}/${data_CASE}/rest/2012-01
if ($do_restarts == true) then
   echo "Restarts starts at "`date`
   cd ${local_arch}/rest
   pwd >>& ${local_arch}/$lists

   # The original ${yr_mo}-* were removed by repack_st_archive 
   # when the tar (into "all types per member") succeeded.
   # The following have been archived to Campaign Storage.
   echo "Restarts ${yr_mo}\*" >>& ${local_arch}/$lists
   ls -d ${yr_mo}* >>& ${local_arch}/$lists
   rm -rf ${yr_mo}*

   # Remove other detritus
   echo "Restarts tar\*.eo" >>& ${local_arch}/$lists
   ls tar*.eo  >>& ${local_arch}/$lists
   rm tar*.eo

   cd ${local_arch}
endif

#------------------------
# Look in 
#    ${data_proj_space}/${data_CASE}/{lnd,atm,ice,rof}/hist/${data_NINST}/${data_CASE}.{clm2,cam,cice,mosart}_*.h[0-9].${year}.nc
# for properly archived files.
# $ ls -l {lnd,atm,ice,rof}/hist/0080/*.{clm2,cam,cice,mosart}_*.h*2012*[cz]
# Or in the following for missed files;
# $ ls {lnd,atm,ice,rof}/hist/*0001.h*00000.nc

if ($do_history == true) then
   set m = 1
   while ($m <= $#components)
      cd ${local_arch}/$components[$m]/hist
      pwd >>& ${local_arch}/$lists
      @ type = 0
      while ($type < 10)
         echo "$components[$m] type $type" >>& ${local_arch}/$lists
         ls ${data_CASE}.$models[$m]_0001.h${type}.${yr_mo}-*.nc > /dev/null
         if ($status != 0) break

         ls ${data_CASE}.$models[$m]_*.h${type}.${yr_mo}-*.nc >>& ${local_arch}/$lists
         rm ${data_CASE}.$models[$m]_*.h${type}.${yr_mo}-*.nc >>& ${local_arch}/$lists
         @ type++
      end
   
      cd ${local_arch}
      @ m++
   end
endif


#------------------------
# Look in ${data_campaign}/${data_CASE}/esp/hist/$yr_mo
# Look in ${data_campaign}/${data_CASE}/atm/hist/$yr_mo
# Look in ${data_campaign}/${data_CASE}/logs/$yr_mo

if ($do_state_space == true) then
   cd ${local_arch}/esp/hist
   pwd >>& ${local_arch}/$lists
   ls -d  $yr_mo/* >>& ${local_arch}/$lists
   rm -rf $yr_mo   >>& ${local_arch}/$lists

   cd ${local_arch}/atm/hist
   pwd >>& ${local_arch}/$lists
   ls -d  $yr_mo/* >>& ${local_arch}/$lists
   rm -rf $yr_mo   >>& ${local_arch}/$lists
   
   # Archive DART log files (and others?)

   cd ${local_arch}/logs
   # This looks misdirected at first, but $lists has 'logs/' in it.
   pwd >>& ${local_arch}/$lists
   ls     $yr_mo >>& ${local_arch}/$lists
   rm -rf $yr_mo >>& ${local_arch}/$lists
   ls     {atm,cpl,ice,lnd,ocn,rof}_00[0-9][02-9].log.*
   rm -rf {atm,cpl,ice,lnd,ocn,rof}_00[0-9][02-9].log.*
   
   cd ${local_arch}

   # Or just the following, but it would give less info if one failed.
#    rm -rf {esp,atm}/hist/$yr_mo logs/$yr_mo  &

endif
#------------------------
if ($do_rundir == true) then
   cd ${local_arch}/../run
   pwd >>& ${local_arch}/lists

   # Remove old inflation restarts that were archived elsewhere
   # and were copied back into rundir by assimilate.csh.
   set files = `ls -t ${data_CASE}.dart.rh.cam_output_priorinf_sd.${yr_mo}*`
   # Skip the most recent.  
   # (The beginning of the next month isn't even in the list.)
   set d = 2
   while ($d <= $#files)
      set date = $files[$d]:r:e
      echo "${data_CASE}.dart.rh.cam*.${date}.*" >>& ${local_arch}/$lists
      ls ${data_CASE}.dart.rh.cam*.${date}.* >>& ${local_arch}/$lists
      rm ${data_CASE}.dart.rh.cam*.${date}.* >>& ${local_arch}/$lists
      @ d++
   end

   # Remove less-than-useful cam_dart_log.${yr_mo}-*.out   
   set files = `ls -t cam_dart_log.${yr_mo}*.out`
   ls $files[2-$#files] >>& ${local_arch}/$lists
   rm $files[2-$#files] >>& ${local_arch}/$lists

   ls finidat_interp_dest_* >& /dev/null
   if ($status == 0) then
      ls finidat_interp_dest_* >>& ${local_arch}/$lists
      rm finidat_interp_dest_* >>& ${local_arch}/$lists
   endif

   cd ${local_arch}
endif
#------------------------

cd archive
gzip $lists

# Wait for all the backrounded 'rm's to finish.
wait

exit

# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

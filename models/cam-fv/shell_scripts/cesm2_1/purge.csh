#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# DART $Id: purge.csh 13192 2019-07-11 21:56:56Z raeder@ucar.edu $

#==========================================================================
# This script removes assimilation output from $DOUT_S_ROOT ($local_archive) after:
# 1) it has been repackaged by repack_st_archive.csh,
# 2) globus has finished moving files to campaign storage.,
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
set do_obs_space   = 'false'
set do_history     = 'true'
set do_state_space = 'true'
set do_rundir      = 'true'

if ($#argv == 0) then
   set project    = /glade/p/nsc/ncis0006/Reanalyses/CAM6_2010
   set campaign   = /gpfs/csfs1/cisl/dares/Reanalyses/CAM6_2010
   set year  = 2011
   set month = 1
   set yr_mo = `printf %4d-%02d ${year} ${month}`
   set obs_space  = Diags_NTrS_${yr_mo}

else if ($#argv == 1) then
   # Request for help; any argument will do.
   echo "Usage:  "
   echo "Before running this script"
   echo "    Run repack_st_archive.csh. "
   echo "    Confirm that it put all the output where it belongs. "
   echo "Call by user or script:"
   echo "   purge.csh project_dir campaign_dir year mo [do_this=false] ... "
   exit

else
   # Script run by another (batch) script or interactively.
   set project  = $1
   set campaign = $2
   set year     = $3
   set month    = $4
   set yr_mo = `printf %4d-%02d ${year} ${month}`
   set obs_space  = Diags_NTrS_${yr_mo}
   # These arguments to turn off parts of the archiving must have the form do_obs_space=false ;
   # no quotes, no spaces.
   if ($?5) set $5
   if ($?6) set $6
   if ($?7) set $7
endif

#--------------------------------------------
set CASEROOT       = $cwd
set CASE           = $CASEROOT:t
set local_arch     = `./xmlquery DOUT_S_ROOT --value`
set ensemble_size  = `./xmlquery NINST_ATM   --value`
set line           = `grep '^[ ]*stages_to_write' input.nml`
set stages_all     = (`echo $line[3-$#line] | sed -e "s#[',]# #g"`)

# "components" = generic pieces of CESM (used in the archive directory names).
# "models" = component instance names (models, used in file names).
set components     = (lnd  atm ice  rof)
set models         = (clm2 cam cice mosart)

#--------------------------------------------
cd $local_arch

set lists = logs/rm_${yr_mo}.lists
if (-f ${lists}.gz) mv ${lists}.gz ${lists}.gz.$$
pwd > $lists

#------------------------
# Look in ${project}/${CASE}/cpl/hist/$INST
if ($do_forcing == true) then
   cd cpl/hist
   pwd >>& ../../$lists

   # The original cpl hist (forcing) files,
   # which have been repackaged into $project.
   foreach type (ha2x3h ha2x1h ha2x1hi ha2x1d hr2x)
      echo "Forcing $type" >>& ../../$lists
      ls ${CASE}.cpl_*.${type}.${yr_mo}-*.nc >>& ../../$lists
      rm ${CASE}.cpl_*.${type}.${yr_mo}-*.nc >>& ../../$lists
   end

   # These files were moved from $project storage to be input to ncrcat.
   # (2020-1-18; that was before implementing Nancy's suggestion to leave
   #  the yearly file on $pr.)
   # An updated version was put into $project, so the locals are not needed.
#    echo "Forcing $year" >>& ../../$lists
#    ls *${year}.nc >>& ../../$lists
#    rm *${year}.nc >>& ../../$lists

   echo "Forcing \*.eo" >>& ../../$lists
   ls *.eo >>& ../../$lists
   rm *.eo >>& ../../$lists

   cd ../..
endif

#------------------------
# Look in ${campaign}/${CASE}/rest/$yr_mo
if ($do_restarts == true) then
   echo "Restarts starts at "`date`
   cd rest
   pwd >>& ../$lists

   # The original ${yr_mo}-* were removed by repack_st_archive 
   # when the tar (into "all types per member") succeeded.
   # The following have been archived to Campaign Storage.
   echo "Restarts ${yr_mo}\*" >>& ../$lists
   ls -d ${yr_mo}* >>& ../$lists
   rm -rf ${yr_mo}*

   # Remove other detritus
   echo "Restarts tar\*.eo" >>& ../$lists
   ls tar*.eo  >>& ../$lists
   rm tar*.eo

   cd ..
endif

#------------------------
# Look in 
#    ${project}/${CASE}/esp/hist/${yr_mo}/${CASE}.dart.cam_obs_seq_final.${yr_mo}.tgz
#    ................................................Diags_NTrS_${yr_mo}.tgz

if ($do_obs_space == true) then
#    cd esp/hist
#    set obs_times_pattern = $yr_mo

   # repack: rm ${CASE}.dart.e.cam_obs_seq_final.${obs_times_pattern}*
   # happened only if the tar file had been moved to $project and had non-0 size.
   # Also; don't do this because the first obs_seq of the next month is in esp/hist,
   #       and we want to keep it.
   # repack; rm ${obs_space} (Diags_NTrS_${yr_mo}) 
   # happened only if the tar file was successfully created and moved to $project.
   #
#    cd ../..
   
endif


#------------------------
# Look in 
#    ${project}/${CASE}/{lnd,atm,ice,rof}/hist/${NINST}/${CASE}.{clm2,cam,cice,mosart}_*.h[0-9].${year}.nc
# for properly archived files.
# $ ls {lnd,atm,ice,rof}/hist/0080/*.{clm2,cam,cice,mosart}_*.h*2010.nc
# Or in the following for missed files;
# $ ls {lnd,atm,ice,rof}/hist/*0001.h*00000.nc

if ($do_history == true) then
   set m = 1
   while ($m <= $#components)
      cd $components[$m]/hist
      pwd >>& ../../$lists
      @ type = 0
      while ($type < 10)
         echo "$components[$m] type $type" >>& ../../$lists
         ls ${CASE}.$models[$m]_0001.h${type}.${yr_mo}-*.nc > /dev/null
         if ($status != 0) break

         ls ${CASE}.$models[$m]_*.h${type}.${yr_mo}-*.nc >>& ../../$lists
         rm ${CASE}.$models[$m]_*.h${type}.${yr_mo}-*.nc >>& ../../$lists
         @ type++
      end
   
      cd ../..
      @ m++
   end
endif


#------------------------
# Look in ${campaign}/${CASE}/esp/hist/$yr_mo
# Look in ${campaign}/${CASE}/atm/hist/$yr_mo
# Look in ${campaign}/${CASE}/logs/$yr_mo

if ($do_state_space == true) then
   cd esp/hist
   pwd >>& ../../$lists
   ls -d  $yr_mo/* >>& ../../$lists
   rm -rf $yr_mo   >>& ../../$lists

   cd ../../atm/hist
   pwd >>& ../../$lists
   ls -d  $yr_mo/* >>& ../../$lists
   rm -rf $yr_mo   >>& ../../$lists
   
   # Archive DART log files (and others?)

   cd ../../logs
   # This looks misdirected at first, but $lists has 'logs/' in it.
   pwd >>& ../$lists
   ls     $yr_mo >>& ../$lists
   rm -rf $yr_mo >>& ../$lists
   ls     {atm,cpl,ice,lnd,ocn,rof}_00[0-9][02-9].log.*
   rm -rf {atm,cpl,ice,lnd,ocn,rof}_00[0-9][02-9].log.*
   
   cd ..

   # Or just the following, but it would give less info if one failed.
#    rm -rf {esp,atm}/hist/$yr_mo logs/$yr_mo  &

endif
#------------------------
if ($do_rundir == true) then
   cd ../run
   pwd >>& ../archive/$lists

   # Remove old inflation restarts that were archived elsewhere
   # and were copied back into rundir by assimilate.csh.
   set files = `ls -t ${CASE}.dart.rh.cam_output_priorinf_sd.${yr_mo}*`
   # Skip the most recent.  
   # (The beginning of the next month isn't even in the list.)
   set d = 2
   while ($d <= $#files)
      set date = $files[$d]:r:e
      echo "${CASE}.dart.rh.cam*.${date}.*" >>& ../archive/$lists
      ls ${CASE}.dart.rh.cam*.${date}.* >>& ../archive/$lists
      rm ${CASE}.dart.rh.cam*.${date}.* >>& ../archive/$lists
      @ d++
   end

   # Remove less-than-useful cam_dart_log.${yr_mo}-*.out   
   set files = `ls -t cam_dart_log.${yr_mo}*.out`
   ls $files[2-$#files] >>& ../archive/$lists
   rm $files[2-$#files] >>& ../archive/$lists

   ls finidat_interp_dest_* >& /dev/null
   if ($status == 0) then
      ls finidat_interp_dest_* >>& ../archive/$lists
      rm finidat_interp_dest_* >>& ../archive/$lists
   endif

   cd ..
endif
#------------------------

cd archive
gzip $lists

# Wait for all the backrounded 'rm's to finish.
wait

exit

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/reanalysis/models/cam-fv/shell_scripts/cesm2_1/purge.csh $
# $Revision: 13192 $
# $Date: 2019-07-11 15:56:56 -0600 (Thu, 11 Jul 2019) $

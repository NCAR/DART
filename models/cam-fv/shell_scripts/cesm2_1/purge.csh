#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#==========================================================================

set MYNAME = $0:t

set INSTRUCTIONS = $MYNAME.$$
cat << ENDOFINSTRUCTIONS >! $INSTRUCTIONS

$MYNAME removes output from \$DOUT_S_ROOT (\$local_archive) for a specified month.
The default action is to simply list what would be purged. An optional argument
triggers the actual removal of the files.

$MYNAME must be run from \$CASEROOT

Usage: 
$MYNAME arg1 [arg2] 

   arg1 may be '--help' , print this help file
   arg1 may be a date in the form YYYY_MM to specify the month to be purged

   If arg2 is specified and has the value 'remove', the files will be removed.

CHECK ALL FINAL ARCHIVE LOCATIONS FOR COMPLETE SETS OF FILES
BEFORE RUNNING $MYNAME

$MYNAME will not allow the current month to be purged.

Before running $MYNAME with the [remove] option ...
 1) Run repack_st_archive.csh
 2) Confirm all the output is where it belongs using pre_purge_check.csh
 3) Run $MYNAME with no arguments to generate a list of files that will
    be removed.
 4) If necessary, edit $MYNAME to change the 
    file types and model/components that will be purged

Example: generate a list of files that would be removed, but do not remove anything
$MYNAME 2014-04

Example: remove files (you will be prompted to continue)
$MYNAME 2014-04 remove

Example: display usage notes
$MYNAME --help

ENDOFINSTRUCTIONS

#==========================================================================

if ($#argv == 1) then
   if ($argv[1] == '--help') then
      cat    $INSTRUCTIONS
      \rm -f $INSTRUCTIONS
      exit 0
   else 
      \rm -f $INSTRUCTIONS
      set purge_date  = `echo $argv[1] | sed -e "s#-# #g"`
      set purge_year  = `echo $purge_date[1] | bc` || exit 1
      set purge_month = `echo $purge_date[2] | bc` || exit 2
      set    SIMPLE_ACTION = 'ls' 
      set RECURSIVE_ACTION = 'ls -R' 
      echo "Performing a 'dry run' ... creating a list of what would be removed."
   endif
else if ($#argv == 2) then
   \rm -f $INSTRUCTIONS
   set purge_date  = `echo $argv[1] | sed -e "s#-# #g"`
   set purge_year  = `echo $purge_date[1] | bc` || exit 1
   set purge_month = `echo $purge_date[2] | bc` || exit 2
   if ($argv[2] == 'remove') then
      set SIMPLE_ACTION = 'rm -v' 
      set RECURSIVE_ACTION = 'rm -rfv' 
      # csh '$<' reads from the keyboard
      echo "WARNING: Removing files for $argv[1]"
      echo "WARNING: Hit any key to continue."
      set go_ahead = $<
   else
      set    SIMPLE_ACTION = 'ls' 
      set RECURSIVE_ACTION = 'ls -R' 
      echo "Performing a 'dry run' ... creating a list of what would be removed."
   endif
else
   cat    $INSTRUCTIONS
   \rm -f $INSTRUCTIONS
   exit 0
endif

#-----------------------------------------------------------------------

if (! -f CaseStatus) then
   echo "ERROR: $MYNAME must be run from the CESM CASEROOT directory"
   exit 1
endif

#-----------------------------------------------------------------------
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

# Get CASE environment variables from the central variables file.
source ./data_scripts.csh
echo "data_CASE  = ${data_CASE}"

#-----------------------------------------------------------------------
# Never purge the current RUNDIR month.

set RUNDIR   = `./xmlquery RUNDIR --value`
set CPL_FILE = `cat ${RUNDIR}/rpointer.drv_0001`
set CPL_DATE = `echo $CPL_FILE:r:e | sed -e "s#-# #g"`
set CPL_YEAR  = `echo $CPL_DATE[1] | bc`
set CPL_MONTH = `echo $CPL_DATE[2] | bc`

set model_yymm = `printf %4d-%02d   ${CPL_YEAR}   ${CPL_MONTH}`
set purge_yymm = `printf %4d-%02d ${purge_year} ${purge_month}`

if ( model_yymm == purge_yymm ) then
   echo "ERROR: cannot purge active month ... $model_yymm"
   echo "ERROR: can only purge previous months."
   exit 2
endif

set yr_mo = $purge_yymm 

#-----------------------------------------------------------------------

set lists_file = ${data_DOUT_S_ROOT}/logs/rm_${yr_mo}.lists
echo "FILE CONTAINING LIST OF FILES TO BE REMOVED:"
echo "   $lists_file"
echo ""

# This will handle running the script in listing mode and then purging mode.
# Both lists_files will be saved.
cd $lists_file:h
if (-f ${lists_file}.gz) mv ${lists_file}.gz ${lists_file}.gz.$$
if (-f ${lists_file}) mv ${lists_file} ${lists_file}.$$

#-----------------------------------------------------------------------
# Record the script settings:

echo "as called: $MYNAME $*"                   >  $lists_file
echo "do_forcing       = ${do_forcing}"        >> $lists_file
echo "do_restarts      = ${do_restarts}"       >> $lists_file
echo "do_history       = ${do_history}"        >> $lists_file
echo "do_state_space   = ${do_state_space}"    >> $lists_file
echo "do_rundir        = ${do_rundir}"         >> $lists_file
echo "components       = ${components}"        >> $lists_file
echo "models           = ${models}"            >> $lists_file
echo "data_CASE        = ${data_CASE}"         >> $lists_file
echo "data_DOUT_S_ROOT = ${data_DOUT_S_ROOT}"  >> $lists_file
echo "SIMPLE_ACTION    = ${SIMPLE_ACTION}"     >> $lists_file
echo "RECURSIVE_ACTION = ${RECURSIVE_ACTION}"  >> $lists_file
echo "running in :"                            >> $lists_file
pwd                                            >> $lists_file

cat $lists_file

#-----------------------------------------------------------------------
# Purge the "forcing" files (cpl history) that came from individual members at single times.
# These have been appended to the yearly files for archiving.
if ($do_forcing == true) then
   cd ${data_DOUT_S_ROOT}/cpl/hist
   pwd >>& $lists_file
   echo "Processing "`pwd`" at "`date`

   # The original cpl hist (forcing) files,
   # which have been repackaged into $data_proj_space.
   foreach type (ha2x3h ha2x1h ha2x1hi ha2x1d hr2x)
      echo "Forcing $type" >>& $lists_file
      ${SIMPLE_ACTION} ${data_CASE}.cpl_*.${type}.${yr_mo}-*.nc >>& $lists_file
   end

   echo "Forcing \*.eo"  >>& $lists_file
   ${SIMPLE_ACTION} *.eo >>& $lists_file

   cd ${data_DOUT_S_ROOT}
endif

#-----------------------------------------------------------------------
# Purge the restart file sets which have been archived to Campaign Storage.
# The original ${yr_mo}-DD-SSSSS directories were removed by repack_st_archive
# when the tar (into "all types per member") succeeded.
# The following ${yr_mo} directories have been archived to Campaign Storage.
if ($do_restarts == true) then
   cd ${data_DOUT_S_ROOT}/rest
   pwd >>& $lists_file
   echo "Processing "`pwd`" at "`date`

   echo      "Restarts ${yr_mo}\*" >>& $lists_file
   ${RECURSIVE_ACTION} ${yr_mo}*   >>& $lists_file

   # Remove other detritus
   echo "Restarts tar\*.eo" >>& $lists_file
   ${SIMPLE_ACTION} tar*.eo >>& $lists_file

   cd ${data_DOUT_S_ROOT}
endif

#-----------------------------------------------------------------------
# Purge component history files (.h[0-9]), 
# which have been tarred into monthly files for each member.
# E.g. {lnd,atm,ice,rof}/hist/0080/*.{clm2,cam,cice,mosart}_*.h*[cz]

if ($do_history == true) then
   set m = 1
   while ($m <= $#components)
      cd ${data_DOUT_S_ROOT}/$components[$m]/hist
      pwd >>& $lists_file
      echo "Processing "`pwd`" at "`date`
      @ type = 0
      while ($type < 10)
         echo "$components[$m] type $type" >>& $lists_file
         ls ${data_CASE}.$models[$m]_0001.h${type}.${yr_mo}-*.nc > /dev/null
         if ($status != 0) break

         ${SIMPLE_ACTION} ${data_CASE}.$models[$m]_*.h${type}.${yr_mo}-*.nc >>& $lists_file
         @ type++
      end
   
      @ m++
   end
endif

#-----------------------------------------------------------------------
# Purge the directories in $DOUT_S_ROOT (scratch ...)
# which have been archived to Campaign Storage.

if ($do_state_space == true) then
   cd ${data_DOUT_S_ROOT}/esp/hist
   echo "Processing "`pwd`" at "`date`
   pwd                          >>& $lists_file
   ${RECURSIVE_ACTION} $yr_mo   >>& $lists_file

   cd ${data_DOUT_S_ROOT}/atm/hist
   echo "Processing "`pwd`" at "`date`
   pwd                                 >>& $lists_file
   ${RECURSIVE_ACTION} $yr_mo          >>& $lists_file
   ${RECURSIVE_ACTION} *.h0.*$yr_mo*   >>& $lists_file
   
   # Archive DART log files (and others?)

   cd ${data_DOUT_S_ROOT}/logs
   echo "Processing "`pwd`" at "`date`
   # This looks misdirected at first, but $lists_file has 'logs/' in it.
   pwd                        >>& $lists_file
   ${RECURSIVE_ACTION} $yr_mo >>& $lists_file
   ${RECURSIVE_ACTION} {atm,cpl,ice,lnd,ocn,rof}_00[0-9][02-9].log.* >>& $lists_file
   
   cd ${data_DOUT_S_ROOT}

endif

#-----------------------------------------------------------------------
# Purge leftover junk in $RUNDIR (scratch ...)
if ($do_rundir == true) then
   cd ${data_DOUT_S_ROOT}/../run
   pwd >>& $lists_file
   echo "Processing "`pwd`" at "`date`

   # Remove old inflation restarts that were archived elsewhere
   # and were copied back into rundir by assimilate.csh.
   set files = `ls -t ${data_CASE}.dart.rh.cam_output_priorinf_sd.${yr_mo}*`
   # Skip the most recent.  
   # (The beginning of the next month isn't even in the list.)
   set d = 2
   while ($d <= $#files)
      set date = $files[$d]:r:e
      ${SIMPLE_ACTION} ${data_CASE}.dart.rh.cam*.${date}.* >>& $lists_file
      @ d++
   end

   # Remove less-than-useful cam_dart_log.${yr_mo}-*.out   
   set files = `ls -t cam_dart_log.${yr_mo}*.out`
   ${SIMPLE_ACTION} $files[2-$#files] >>& $lists_file

   ls finidat_interp_dest_* >& /dev/null
   if ($status == 0) then
      ${SIMPLE_ACTION} finidat_interp_dest_* >>& $lists_file
   endif

   cd ${data_DOUT_S_ROOT}
endif

#-----------------------------------------------------------------------

cd ${data_DOUT_S_ROOT}/logs
gzip $lists_file

exit 0

#!/bin/csh -f
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# This script is used to archive DART and CESM files when running multiple
# CESM/DART cycles in a single batch job.  See the README file for more details.

#BSUB -o archive_cycles.%J
#BSUB -e archive_cycles.%J
#BSUB -W 11:00
#BSUB -N 
#BSUB -q caldera
#BSUB -n 1
#BSUB -P ZZZZZZZZ
#BSUB -J archive_cycles.csh

# ==============================================================================
# Load environment variables from CESM 
# ==============================================================================

cd BOGUS_CASE

source ./Tools/ccsm_getenv || exit -1

# Set the archive frequency for restart sets. The others will be removed.
set archive_Nth_days = 4

# ==============================================================================
# standard commands:
# 
# If you are running on a machine where the standard commands are not in the
# expected location, add a case for them below.
# ==============================================================================

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set REMOVE = '/usr/local/bin/rm -fr'
      set   LINK = '/usr/local/bin/ln -fvs'
      breaksw

   default:
      # NERSC "hopper", NWSC "yellowstone"
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set REMOVE = '/bin/rm -fr'
      set   LINK = '/bin/ln -fvs'
      breaksw
endsw

# ==============================================================================
# Use the CESM coupler files to make a file list; options are to
#  delete them, archive them, or leave them in place. 
# ==============================================================================

cd $RUNDIR
set most_recent = `ls -1t *cpl*.r* | head -1`
set most_recent_date =  `echo $most_recent | sed "s/\.nc//; s/^.*\.r\.//;"`
echo THE MOST RECENT DATE, $most_recent_date, RESTARTS WILL be saved 

foreach cplfile ( *cpl*.r* )

   # Extract the date from the cpl restart file name
   set date = `echo $cplfile | sed "s/\.nc//; s/^.*\.r\.//;"`
   
   if ($date != $most_recent_date) then
      set temp = `echo $date | sed -e "s#-# #g"`
      set day = $temp[3]
      set secs = $temp[4]
      if ( $secs == 00000 && (`expr $day % $archive_Nth_days` == 0 || $day == 28) ) then

          echo $date : archiving files from this date
          if ( -e $DOUT_S_ROOT/rest/${date} ) then
             echo directory already exists
          else
             echo creating directory
             mkdir -p $DOUT_S_ROOT/rest/${date}
          endif
          $MOVE *.r*.${date}*                  $DOUT_S_ROOT/rest/${date}
          $MOVE *.i.${date}*                   $DOUT_S_ROOT/rest/${date}
          $MOVE *prior_inflate_restart.${date} $DOUT_S_ROOT/rest/${date}
          $MOVE *post_inflate_restart.${date}  $DOUT_S_ROOT/rest/${date}
          $COPY *cam*.h*.${date}*              $DOUT_S_ROOT/rest/${date}
          $COPY *clm*.h*.${date}*              $DOUT_S_ROOT/rest/${date}
      else
         echo $date : deleting files from this date
         $REMOVE *.r*.${date}* 
         $REMOVE *.i.${date}*  
         $REMOVE *prior_inflate_restart.${date}
         $REMOVE *post_inflate_restart.${date}
      endif
   
   else
       echo $date : preserving files from this date
   endif

end

# ==============================================================================
# archive history files
# ==============================================================================

# Save most recent day's worth of history files, for potential continuing accumulation.
# This assumes you are doing 6 hour assimilation cycles.
set times = (00000 21600 43200 64800)
foreach t ($times)
   set most_recent = `ls -1t *cpl*.r*-$t* | head -1`
   if ($status != 0) then
      # Move on to the next iteration of the t loop.
      continue
   endif

   set most_recent_date =  `echo $most_recent | sed "s/\.nc//; s/^.*\.r\.//;"`
   echo THE MOST RECENT DATE, $most_recent_date, HISTORY FILES WILL be saved 
   
   if ( -e TEMP ) then
      echo TEMP directory already exists to hold current history files
   else
      echo creating TEMP directory to hold current history files
      mkdir TEMP 
   endif

   #move current history files into this directory
   $MOVE *cam*.h*.${most_recent_date}*   TEMP
   $MOVE *pop*.h*.${most_recent_date}*   TEMP
   $MOVE *pop*.d*.${most_recent_date}*   TEMP
   $MOVE *rtm*.h*.${most_recent_date}*   TEMP
   $MOVE *clm2*.h*.${most_recent_date}*  TEMP
   $MOVE *cice*.h*.${most_recent_date}*  TEMP
   # Move CAM-SE grid files out of the way
   $MOVE *Mapping*.nc                    TEMP
end

# Now archive all other history files.
# All times (except for those hidden in TEMP) will be moved to the archive directory.
# Excluding history files based on times requires extra code here.

# Each of these entries will have * prepended and appended to it.
set files = ('cpl.log.'       'cesm.log.'      'dart_log.'      'P*Diag'       'True'            \
             'obs_seq'        '$CASE.cam*.h'   'atm_0001*.log.' '$CASE.clm*.h' 'lnd_0001*.log.'  \
             'ice_0001*.log.' 'atm.log.'       'lnd.log.')

#  These are parallel lists; entries here must correspond exactly to the file list directly above.
set dests = ('cpl/logs'       'cpl/logs'       'dart/logs'      'dart/hist'      'dart/hist' \
             'dart/hist'      'atm/hist'       'atm/logs'       'lnd/hist'       'lnd/logs'  \
             'ice/logs'       'atm/logs'       'lnd/logs')


if ($#files != $#dests) then
   echo "Wordlists 'files' and 'dests' must have the same number of words in them"
   exit 89
endif

# Make copies of the obs_seq.final files in an unarchived place,
# to make obs space diagnostics easier.
set o_s_finals = $RUNDIR:h/Obs_seqs
if (! -d ${o_s_finals}) mkdir ${o_s_finals}

set f = 1
while ($f <= $#files)
   set file_set = "*$files[$f]*"
   ls $file_set >& /dev/null 
   if ($status != 0) then
      echo "Finished with all files matching $files[$f]"
   else
      if ("$files[$f]" == 'obs_seq' ) then
         $COPY $file_set ${o_s_finals}
      endif
      
      if (! -d $DOUT_S_ROOT/$dests[$f] ) then
         echo "Making $DOUT_S_ROOT/$dests[$f] " 
         mkdir -p $DOUT_S_ROOT/$dests[$f]
      endif

      $MOVE $file_set $DOUT_S_ROOT/$dests[$f]
   endif 
   @ f++
end

# Remove all log files that haven't been archived.

$REMOVE *log*


# move the history files in TEMP back into RUNDIR:
$MOVE TEMP/* .

# ==============================================================================
# run the long term archiver if requested
# ==============================================================================

if ($DOUT_L_MS == 'TRUE') then
   cd $DOUT_S_ROOT
   $CASEROOT/Tools/lt_archive.sh -m copy_dirs_hsi
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


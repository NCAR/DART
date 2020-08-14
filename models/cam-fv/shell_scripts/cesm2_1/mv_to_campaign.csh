#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download


# Script to package output from a CAM+DART assimilation,
# especially the Reanalysis project (2019),
# and send it to Campaign storage.
# It's derived from ./mv_to_campaign.sample.csh.
# Documentation of that script is in 
#   https://www2.cisl.ucar.edu/sites/default/files/CISL_GlobusCLI_Nov2018.html

# To get access to the python globus command it seems necessary 
# to do the following from the command line (not in python):

# > module load gnu python    (even though the python command is available without loading)
# > ncar_pylib
# > globus --help

# This will show the commands that can be given to globus, like 'endpoint'
# The subcommands of 'endpoint', like 'search' can be seen with

# > globus list-commands

# But that shows subcommands of all commands.  Focus the search with 

# > globus [command [subcommand]] --help
# For example 
# > globus endpoint activate --help

# will show obscura such as --myproxy (below).

# >> Add error checking on arguments
# and a "Usage:" section

# Case name
set CASENAME = "$1"
# or set these explicitly in here?
# Depends whether it will be run from another script.

# Set analysis time.
# But we'll be archiving multiple times at once.
# But maybe by month, so YYYY-MM would be useful.
set TIMESTR="$2"

# Need log file time stamp time?
set STAMP = "$3"

# Beware: the label cannot contain slashes.  valid chars are only:
# A-Z, a-z, 0-9, space, hyphen, underscore, and comma
# No '.'?
# Max length is 128 chars.
set LABEL = "copy of $CASENAME dares project files for $TIMESTR" 

# Declare paths to use in script (EDIT THESE BEFORE RUNNING!)
# set SRC_DIR=/glade/scratch/${USER}/${CASENAME}/archive
set SRC_DIR=/glade/p/cisl/dares/Reanalyses/${CASENAME}/archive
# Campaign Storage
# Change a dir name to Reanalyses or ...sis?
set CS_DIR=/gpfs/csfs1/cisl/dares/Reanalysis/CAM6_2017/

# Done with input parameters.
#=================================================================

set AnY=`date -d $TIMESTR '+%Y'`
set AnM=`date -d $TIMESTR '+%m'`
# set AnD=`date -d $TIMESTR '+%d'`
set AN_DATE = ${AnY}${AnM}
set glog   = globus_${AN_DATE}.log

# start with an empty log
cd $SRC_DIR:h
rm -f $glog globus-batch-dirs.txt globus-batch-files.txt
echo Copy $SRC_DIR to campaign storage $CS_DIR >>& $glog

# Load Python to get the CLI
module load gnu python
# Activate the NCAR Python Library (NPL) virtual environment 
# for version given as argument.
# This command activates the 'globus' command, used below.
ncar_pylib 20190118

# Retrieve endpoint IDs and store them as variables
# Access to the globus command comes through the python module.
# That module requires  ncarenv/1.2  gnu/6.3.0  ncarcompilers/0.4.1  
# gnu replaces intel/#.#.# that I have already loaded.
# OK because that load expires with the end of this script.
# --filter-owner-id not documented.
# --jq              is A JMESPath expression to apply to json output.
#                   Takes precedence over any specified '--format'
# But this has a '--format UNIX' anyway.
# 
# EP = endpoint
set EP_SRC=`globus endpoint search 'NCAR GLADE'            \
                --filter-owner-id ncar@globusid.org         \
                --jq 'DATA[0].id' --format UNIX`
set EP_CS=`globus endpoint search 'NCAR Campaign Storage' \
                --filter-owner-id ncar@globusid.org         \
                --jq 'DATA[0].id' --format UNIX`

# Nancy had to add this activation before being able to see 
# the EP_CS location using globus.
# (E.g. > globus ls ${EP_CS}:/gpfs/csfs1/cisl/dares/Reanalysis/).
# It seems to activate the endpoints without requiring a password.
# (Or was she already logged into globus from previous commands?)
foreach ep ($EP_SRC $EP_CS)
   # Check if endpoint is activated
   # (we don't care about output, only return code)
   globus endpoint is-activated $ep >& /dev/null
   if ( $status != 0 ) then
      globus endpoint activate --myproxy --myproxy-lifetime 1 $ep
      if ( $status != 0 ) then
         echo "Fatal: NCAR endpoint $ep isn't activated." > $glog
         echo "Aborting transfer..." >> $glog
         echo "Failed: $AN_DATE to Campaign Storage!" > ~/GLOBUS-ERROR.$AN_DATE
         exit 1
      endif
   endif
end

set EXPIRE=`globus endpoint is-activated                \
            --jq expire_time -F unix $EP_SRC`
echo "NCAR endpoints active until $EXPIRE" > $glog

# Check if destination directory exists; if not, create it
globus ls ${EP_CS}:$CS_DIR >& /dev/null

if ( $status != 0 ) then
    globus mkdir ${EP_CS}:$CS_DIR >>& $glog
else
    echo $CS_DIR already exists on campaign store >>& $glog
endif

set DESTDIR=${CS_DIR}/$AN_DATE
globus mkdir ${EP_CS}:${DESTDIR} >>& $glog

# Create a list of files to archive.
# This is mysterious, since I don't know how the original file names look.
#    set BATCHAnMT="${RUNDIR}/\1 ${DESTDIR}/\1"
#    ls -1 fcst*.nc | sed "s|\(.*\)|${BATCHAnMT}|" > globus-batch.txt
# Hopefully the output to globus-batch.txt does not have special formatting.
# It may:
# "The batch file needs to have full source and destination paths
#  We use a sed command to format our ls output and store it as a bash variable SOUT"
# 
# The list will depend on the packaging of files:
# laptop:/Users/raeder/DAI/ATM_forcXX/CAM6_setup/campaign_storage

# this finds directories as well as files.  not sure what we need here.
#ls -1R . | sed -e "s;.*;${SRC_DIR}/& ${CS_DIR}/&;" > globus-batch.txt
# find . -type d | sed -e "s;.*;globus mkdir ${EP_CS}:${CS_DIR}/&;" >! globus-batch-dirs.txt

# possibly this:

# cd into the source dir
cd $SRC_DIR

# find all dirs, and get rid of the ./ at the start of each subdir name
find . -type d | sed -e "s;..;;"  >! ../globus-batch-dirs.txt

# find all files, get rid of the ./ at the start of each filename, and 
# convert them to 2 full pathnames: the source and the destination
find . -type f | sed -e "s;..\(.*\);${SRC_DIR}/\1 ${DESTDIR}/\1;" >! ../globus-batch-files.txt

# do we need to check for their existance first? doing so to be safe.
echo Creating needed subdirectories
foreach SUBDIR ( `cat ../globus-batch-dirs.txt` )
   set target = ${EP_CS}:${CS_DIR}/${SUBDIR}

   # Check if destination directory already exists; if not, create it
   globus ls $target >& /dev/null

   if ( $status != 0 ) then
      echo Making $target on campaign store   >>& ../$glog
      globus mkdir $target                    >>& ../$glog
   else
      echo Subdir $target already exists on campaign store >>& ../$glog
   endif
end

echo files to be copied are in globus-batch-files.txt

# Start copy of GLADE data holdings to CS
# Finally, we use the variable contents as stdin to our globus batch transfer
globus transfer $EP_SRC $EP_CS                           \
    --label "$LABEL"              \
    --batch < ../globus-batch-files.txt >>& ../$glog

echo ""
echo Output of this script is in $SRC_PARENT/$glog.
echo Transfer is asynchronous.  If successfully started, 
echo you will receive email when it is complete

echo ""
echo Ending script to copy the contents of $SRC_DIR to campaign storage at `date`

exit 0



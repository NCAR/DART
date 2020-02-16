#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# $Id$

# Script to send output from a CAM+DART assimilation,
# especially the Reanalysis project (2019), to Campaign Storage.
# >>> Before running this script:
#     1) (re)package the files into a directory containing only files to be archived.
#        The final subdirectory is often a CESM style date string; YYYY-MM-DD-SSSSS
#     2) You may need to interactively log in to globus.
#        See https://www2.cisl.ucar.edu/resources/storage-and-file-systems/globus-file-transfers
#        for details.
#        But may need to first
#        > ssh data-access.ucar.edu 
#          (alternative: ssh username@data-access.ucar.edu)
#          or open a window on casper.
#        It may be enough to 
#        > module load gnu python
#        > ncar_pylib 20190326
#        OR log in through the globus file manager.
#        In either case:
#        > globus login
#        This seems unnecessary:
#           Copy the resulting URL into your browser
#           Log in there to get the access code.
#           Enter the code into the interactive globus prompt.

# This script was derived from ./mv_to_campaign.sample.csh.
# Documentation of that script is in 
#   https://www2.cisl.ucar.edu/sites/default/files/CISL_GlobusCLI_Nov2018.html

if ($#argv == 1) then

   # Request for help; any argument will do.
   echo "Usage: call by user or script:"
   echo "   mv_to_campaign.csh CASE TIME_STR SRC_DIR CS_DIR"
   echo '      CASE     = CESM $CASE name.'
   echo "      TIME_STR = CESM format time string associated with the data; YYYY-MM-DD-SSSSS"
   echo "      SRC_DIR  = the directory to archive in Campaign Storage"
   echo "      CS_DIR   = Campaign Storage"
   exit

else if ($#argv == 0) then

   set CASE     = Test_mv_to_globus
   # set CASE     = f.e21.FHIST_BGC.f09_025.CAM6assim.003
   set TIME_STR = 2019-04-23-91800
   # set TIME_STR = 2017-01-02-00000
   set SRC_DIR  = /glade/scratch/${USER}/${CASE}/${TIME_STR}
   # set SRC_DIR  = /glade/scratch/${USER}/${CASE}/archive/rest/${TIME_STR}
   # set SRC_DIR  = /glade/p/cisl/dares/Reanalyses/CAM6_2017/${CASE}/archive
   set CS_DIR   = /gpfs/csfs1/cisl/dares/Reanalyses/CAM6_2017/${CASE}

else if ($#argv == 4) then

   set CASE     = $1 
   set TIME_STR = $2 
   set SRC_DIR  = $3
   set CS_DIR   = $4

else
   
   echo "ERROR: This script requires 0, 1, or 4 arguments"
   exit
endif

# Done with input parameters.
#=================================================================

# Beware: the label cannot contain slashes.  valid chars are only:
# A-Z, a-z, 0-9, space, hyphen, underscore, and comma
# No '.'; replace with ','
# Max length is 128 chars.
set clean_CASE = `echo $CASE | sed -e "s#\.#,#g"`  
set LABEL = "copy of $clean_CASE dares project files for $TIME_STR" 

set AN_DATE = $SRC_DIR:t

cd $SRC_DIR:h

# start with an empty log
set glog = globus_${AN_DATE}_$$.log
rm -f $glog globus-batch-dirs.txt globus-batch-files.txt
echo Copy $SRC_DIR to campaign storage $CS_DIR >>& $glog

# Load Python to get the globus Command Line Interface.
# Learn more about globus using
# > globus --help
# This will show the commands that can be given to globus, like 'endpoint'
# The subcommands of 'endpoint', like 'search' can be seen with
# > globus list-commands
# But that shows subcommands of all commands.  Focus the search with 
# > globus [command [subcommand]] --help
# For example 
# > globus endpoint activate --help
# will show obscura such as --myproxy (below).
module load gnu python

# Activate the NCAR Python Library (NPL) virtual environment 
# for the version given as the argument (no arg = use default
# which doesn't work 2019-6-21).
# This command activates the 'globus' command, used below.
ncar_pylib 20190326
# ncar_pylib 

# Retrieve endpoint IDs and store them as variables using globus.
# --filter-owner-id not documented.
# --jq              is A JMESPath expression to apply to json output.
#                   Takes precedence over any specified '--format'
# But this has a '--format UNIX' anyway.
# 
# EP = endpoint
set EP_SRC = `globus endpoint search 'NCAR GLADE'           \
                --filter-owner-id ncar@globusid.org         \
                --jq 'DATA[0].id' --format UNIX`
set EP_CS = `globus endpoint search 'NCAR Campaign Storage' \
                --filter-owner-id ncar@globusid.org         \
                --jq 'DATA[0].id' --format UNIX`

echo EP_SRC = $EP_SRC
echo EP_CS = $EP_CS

# Add these activations in order to use (and see) the end points.
# (E.g. > globus ls ${EP_CS}:/gpfs/csfs1/cisl/dares/Reanalysis/).
# It seems to activate the endpoints without requiring a password.
# That's because the credentials were set up manually 
# with a lifetime of 30 days.  "activate" first tries auto-activation.
# If that fails, it tries the activation method from the command line.
# But all of those options require a browser (not available on batch nodes)
# or user name and password (which is insecure if provided by the script).
# So this activation relies on auto-activation.
# It may be redundant, but maybe this script needs an explicit activation.
foreach ep ($EP_SRC $EP_CS)
   # Check if endpoint is activated
   # (we don't care about output, only return code)
   # globus endpoint is-activated $ep >& /dev/null
   set EXPIRE = `globus endpoint is-activated \
                 --jq expire_time -F unix $ep`
   # if ( $status != 0 ) then
   if ( `echo $EXPIRE` == 'None' ) then
      globus endpoint activate $ep
      if ( $status != 0 ) then
         echo "Fatal: NCAR endpoint $ep isn't activated." > $glog
         echo "       Aborting transfer..." >> $glog
         echo "Failed: $AN_DATE to Campaign Storage!" > ~/GLOBUS-ERROR.$AN_DATE
         exit 1
      endif
   else
      echo Endpoint $ep is activated until $EXPIRE
      echo "NCAR endpoint $ep active until $EXPIRE" > $glog
   endif
end

set DEST_DIR = ${CS_DIR}/$AN_DATE

# Cd from where task files will be created into the source dir where the data files are.
cd $SRC_DIR

# Start copy of GLADE data holdings to Campaign Storage.
# --recursive  to duplicate an entire directory and all its contents. 
# --sync-level mtime anything with a newer modification time gets copied,
#              If a transfer fails, CHECKSUM must be used to restart the transfer. 
#              > > > All other levels can lead to data corruption. < < <
# ?            Is restarting different than submitting another transfer request?
# ?            If a transfer partially succeeds, does the partial file have a newer
#              date than the old, so that mtime will prevent the new transfer?
# --batch < ../globus-batch-files.txt   NOT needed when moving whole directories\
globus transfer                         \
    --recursive --sync-level mtime      \
    --label "$LABEL"                    \
    ${EP_SRC}:${SRC_DIR} ${EP_CS}:${DEST_DIR}   >>& ../$glog


echo ""
echo "Transfer is asynchronous."
echo "IF successfully started, you will receive email when it is complete."
echo "CHECK $SRC_DIR:h/$glog."

echo ""
echo "Ending script to copy the contents of $SRC_DIR "
echo "to campaign storage at `date`"

exit 0


# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

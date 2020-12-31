#!/bin/csh

# Run this after gen_rean_diags creates the pictures.
# -a, --archive               archive mode; 
# -v, --verbose               increase verbosity (default is -q quiet)
# -z, --compress              compress file data during the transfer
# -C, --cvs-exclude           auto-ignore files in the same way CVS does
# -e, --rsh=COMMAND           specify the remote shell to use
#     --exclude=PATTERN       exclude files matching PATTERN

set my_pth = `pwd`
set diags_dir = $my_pth:t

# Archive the ps, pdf, and obs_diag_output.nc files 
# while the png files are sent to the web site.
tar -z -c -f ../${diags_dir}.tgz [^w]*  &

# The trailing slashes are important in this syntax

set MACHINE=webpub.ucar.edu
set MACHINEDIR=/test/image/pub/DART/Reanalysis
# The web page where goes first is
# https://test.www.image.ucar.edu/pub/DART/Reanalysis
# but I must be connected via VPN.

# The last directory in LOCALDIR must be the one created by "gen_rean_diags.m" 
# Make sure there is no trailing '/' on LOCALDIR or MACHINEDIR
# Orig:
# set LOCALDIR=~/Documents/NCAR/Reanalysis/Diags_Feb-2011/web_Feb-2011
set web_dir = `ls -d web_*`
if ($status != 0) then
   echo "ERROR: no web_ directory"
   exit
endif
set LOCALDIR = ${my_pth}/${web_dir}

# ? Exclude git files instead of .svn?
echo "rsync -avz --rsh ssh --exclude '"'.svn*'"'\\" 
echo "${LOCALDIR} ${MACHINE}:${MACHINEDIR}"

rsync -avz --rsh ssh --exclude '.svn*' ${LOCALDIR} ${MACHINE}:${MACHINEDIR}

wait

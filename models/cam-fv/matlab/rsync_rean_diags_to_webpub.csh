#!/bin/csh

# -a, --archive               archive mode; 
# -v, --verbose               increase verbosity
# -z, --compress              compress file data during the transfer
# -C, --cvs-exclude           auto-ignore files in the same way CVS does
# -e, --rsh=COMMAND           specify the remote shell to use
#     --exclude=PATTERN       exclude files matching PATTERN

# The trailing slashes are important in this syntax

set MACHINE=webpub.ucar.edu
set MACHINEDIR=/test/image/pub/DART/Reanalysis

# The last directory in LOCALDIR must be the one created by "gen_rean_diags.m" 
# Make sure there is no trailing '/' on LOCALDIR or MACHINEDIR
set LOCALDIR=~/Documents/NCAR/Reanalysis/Diags_Feb-2011/web_Feb-2011

rsync -avz --rsh ssh --exclude '.svn*' ${LOCALDIR} ${MACHINE}:${MACHINEDIR}


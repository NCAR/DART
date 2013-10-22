#!/bin/csh 
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# run obs_diag on the mainframe a month at a time.  produce
# obs_diag_output.nc file(s) which need to be copied to 
# another machine which has matlab to finish generating 
# the plots that are normally looked at.

# this program is expected to be started from the directory 
# which has a collection of obs_NNNN dirs, each of which 
# contains an obs_seq.final file.

#-----------------------------------------------------------------------------
# run in share queue on bluefire.
#-----------------------------------------------------------------------------

#BSUB -J obs_diags
#BSUB -o obs_diags.%J.log
#BSUB -P xxxxxxxx
#BSUB -q share
#BSUB -W 2:00
#BSUB -n 1


# general settings - check these to see they are appropriate

# where to find the DART obs_diag program, and a template for
# the input.nml file.   obs_diag should exist in this dir.
set DIAG = $DART/models/cam/work


# time period our standard experiments have covered.  if you
# move to another time, you'll need to regenerate an appropriate
# table.  you can set the start and stop month number below.
# remember month 1 is the spinup month, so generally our
# experiments run from month 2 to 13.

set startmon = 2
set endmon   = 3


# tables intended to ease the calendar computations.
# these are very date specific and MUST be changed for any other
# time spans that are run.  day 1 = aug 1, 2006
# should go one month past last day of simulation
set months = (aug06 sep06 oct06 nov06 dec06 jan07 feb07 mar07 \
              apr07 may07 jun07 jul07 aug07 sep07)

#               1  2  3   4   5   6   7   8   9  10  11  12  13  14
set startd  = ( 1 32 62  93 123 154 185 213 244 274 305 335 366 397)
set endd    = (31 61 92 122 153 184 212 243 273 304 334 365 396 426)
set year    = (06 06 06  06  06  07  07  07  07  07  07  07  07  07)
set calmon  = ( 8  9 10  11  12   1   2   3   4   5   6   7   8   9)
set daysmon = (31 30 31  30  31  31  28  31  30  31  30  31  31  30)

@ monnum = $startmon
while ( $monnum <= $endmon )
  set diagname = $months[$monnum]

  # assume there is an input.nml.diag.template here with the right 
  # patterns for sed to search and replace below.
  # should error check for this file here.

  set sdy = `printf '%04d' $startd[$monnum]`
  set smn = $calmon[$monnum]
  set mdy = $daysmon[$monnum]
  set syr = $year[$monnum]
   
  sed -e "s/NNNN/$sdy/" \
      -e "s/YY/$syr/"    \
      -e "s/MM/$smn/"    \
      -e "s/DD/$mdy/"    \
      input.nml.diag.template >! input.nml

  echo updated input.nml here

  echo running obs_diag here
  $DIAG/obs_diag >&! obs_diag.out

  set nextdir = ${startd[$monnum]}-${endd[$monnum]}s0_def_reg

  mkdir -p $nextdir

  touch $diagname
  mv -f $diagname obs_diag_output.nc input.nml obs_diag.out \
        LargeInnov.txt $nextdir

  @ monnum++
end

echo all done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


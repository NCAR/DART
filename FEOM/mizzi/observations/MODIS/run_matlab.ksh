#!/bin/ksh -x
#
# Script to run matlab
bsub -Is -q geyser -W 4:00 -n 1 -P NACD0002 matlab
exit

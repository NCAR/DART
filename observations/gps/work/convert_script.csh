#!/bin/csh
#
# Main script:
# generate multiple days of gps observations
#
# calls the cosmic_to_obsseq script with a date, the working
# directory location, and whether to download the data automatically
# from the cosmic web site (downloading data requires signing up 
# for a username/password to access the site, and then setting 
# the username and password here before running this script.)
#

setenv cosmic_user xxx
setenv cosmic_pw   yyy

./cosmic_to_obsseq.csh 20061101 .. yes
./cosmic_to_obsseq.csh 20061102 .. yes
./cosmic_to_obsseq.csh 20061103 .. yes
./cosmic_to_obsseq.csh 20061104 .. yes
./cosmic_to_obsseq.csh 20061105 .. yes
./cosmic_to_obsseq.csh 20061106 .. yes
./cosmic_to_obsseq.csh 20061107 .. yes

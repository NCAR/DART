#!/bin/csh
#
# Main script:
# generate multiple days of gps observations
#
# calls the cosmic_to_obsseq script with 4 args:
#
#  - the date in YYYYMMDD format
#  - the working directory location
#  - whether to download the data automatically from the cosmic web site 
#     (downloading data requires signing up for a username/password to 
#     access the site, and then setting the username and password here 
#     before running this script.)  set to no if the data has already been
#     downloaded separately before now.
#  - whether to delete the data automatically from the local disk after the
#     conversion is done.
#

setenv cosmic_user xxx
setenv cosmic_pw   yyy

./cosmic_to_obsseq.csh 20061101 ../cosmic no no
./cosmic_to_obsseq.csh 20061102 ../cosmic no no
./cosmic_to_obsseq.csh 20061103 ../cosmic no no
./cosmic_to_obsseq.csh 20061104 ../cosmic no no
./cosmic_to_obsseq.csh 20061105 ../cosmic no no
./cosmic_to_obsseq.csh 20061106 ../cosmic no no
./cosmic_to_obsseq.csh 20061107 ../cosmic no no

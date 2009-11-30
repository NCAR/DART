#!/bin/csh
#
# Main script:
# download multiple days of gps files
#
# calls the cosmic_download script with a date and the working
# directory location.  downloading data requires signing up 
# for a username/password to access the site, and then setting 
# the username and password here before running this script.
#

#BSUB -J get_gps
#BSUB -o gps.%J.log
#BSUB -q standby
#BSUB -n 1
#BSUB -W 12:00

setenv cosmic_user xxx
setenv cosmic_pw   xxx

./cosmic_download.csh 20071001 ../cosmic
./cosmic_download.csh 20071002 ../cosmic
./cosmic_download.csh 20071003 ../cosmic
./cosmic_download.csh 20071004 ../cosmic
./cosmic_download.csh 20071005 ../cosmic
./cosmic_download.csh 20071006 ../cosmic
./cosmic_download.csh 20071007 ../cosmic
 
./cosmic_download.csh 20071008 ../cosmic
./cosmic_download.csh 20071009 ../cosmic
./cosmic_download.csh 20071010 ../cosmic
./cosmic_download.csh 20071011 ../cosmic
./cosmic_download.csh 20071012 ../cosmic
./cosmic_download.csh 20071013 ../cosmic
./cosmic_download.csh 20071014 ../cosmic
 
./cosmic_download.csh 20071015 ../cosmic
./cosmic_download.csh 20071016 ../cosmic
./cosmic_download.csh 20071017 ../cosmic
./cosmic_download.csh 20071018 ../cosmic
./cosmic_download.csh 20071019 ../cosmic
./cosmic_download.csh 20071020 ../cosmic
./cosmic_download.csh 20071021 ../cosmic
 
./cosmic_download.csh 20071022 ../cosmic
./cosmic_download.csh 20071023 ../cosmic
./cosmic_download.csh 20071024 ../cosmic
./cosmic_download.csh 20071025 ../cosmic
./cosmic_download.csh 20071026 ../cosmic
./cosmic_download.csh 20071027 ../cosmic
./cosmic_download.csh 20071028 ../cosmic

./cosmic_download.csh 20071029 ../cosmic
./cosmic_download.csh 20071030 ../cosmic
./cosmic_download.csh 20071031 ../cosmic

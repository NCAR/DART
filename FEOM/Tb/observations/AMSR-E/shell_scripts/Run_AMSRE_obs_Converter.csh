#!/bin/tcsh -f
#
# Download the AMSRE data with:
#
# wget -t 0 -T 120 -c -r -nH --cut-dirs=4 \
#    ftp://sidads.colorado.edu/pub/DATASETS/nsidc0301_amsre_ease_grid_tbs/global/2003/

# Directories for Ally
# set RAW_OBS_DIR = /gpfsm/dnb52/projects/p31/DATA/DART_DATA/Observations/AMSRE_obs/AMSRE_RAW
# set AMSRE_WORKDIR = /gpfsm/dnb31/atoure/DART/observations/AMSR-E

# Directories for Tim
set   RAW_OBS_DIR = /Users/thoar/Desktop/EASE_Grid/2011_north
set AMSRE_WORKDIR = /Users/thoar/Desktop/EASE_Grid/temp
set      DART_DIR = /Users/thoar/svn/DART/Tb/observations/AMSR-E/work

set nonomatch

# If the work directory does not exist, make it.
if  ( ! -d ${AMSRE_WORKDIR} ) then
  mkdir -p ${AMSRE_WORKDIR}
endif

cd ${AMSRE_WORKDIR}

\cp -v ${DART_DIR}/input.nml        . || exit 1
\cp -v ${DART_DIR}/ease_grid_to_obs . || exit 2
\cp -v ${DART_DIR}/advance_time     . || exit 3

# Set DOY equal to the first day-of-year of interest.
# Set NDAYS equal to the number of days to convert.
# Keep in mind it takes a LONG time to convert 1 day of
# all frequencies, polarizations, passes. For the NH alone,
# this takes more than 40 minutes and results in 6.3+ MILLION
# observations - for a single day! 

@ YEAR = 2011
@ DOY = 5
@ NDAYS = 4
set FREQ = 89V
set FBASE = ID2r3-AMSRE-NL

# do not change anything below this line

@ IDAY = 1

while ($IDAY <= $NDAYS)

   set FileBaseDate = `printf %s%04d%03d ${FBASE} $YEAR $DOY`

   echo "looking for ${RAW_OBS_DIR}/${YEAR}/${FileBaseDate}*"
   \ls -1 ${RAW_OBS_DIR}/${YEAR}/${FileBaseDate}* >! file_list.$$

   # subset just the frequency of interest 

   grep $FREQ file_list.$$ >! file_list.txt

   if ( ! -z file_list.txt ) then

       echo "Input files are:"
       cat file_list.txt

       # The date format of the output file name has to be yyyy-mm-dd-00000
       # (in part, because the files always pertain to midnight)
       # We must convert the day-of-year to month-day since the files are 
       # stored in directories named YYYYMM

       set FILE_YYYYMMDD = `echo ${YEAR}010100 +${DOY}d-1d -f ccyy-mm-dd | ./advance_time`
       set FILE_YYYYMM   = `echo ${YEAR}010100 +${DOY}d-1d -f ccyy-mm    | ./advance_time`
       set YYYYMM = `echo $FILE_YYYYMM | sed 's/-//'`

       set Output_Obs_dir = "DART_OBS_SEQ/${YYYYMM}"
       echo "The processed data dir is " $Output_Obs_dir

       if  ( ! -d ${RAW_OBS_DIR}/${Output_Obs_dir} ) then
         mkdir -p ${RAW_OBS_DIR}/${Output_Obs_dir}
       endif

       # specify output Output_fileName, format: obs_seq.yyyy-mm-dd-00000.out
       set Date_output_file = ${FILE_YYYYMMDD}-00000
       set Output_fileName = "obs_seq.${Date_output_file}.out"
       set Full_Output_fileName = "${RAW_OBS_DIR}/${Output_Obs_dir}/${Output_fileName}"
       echo "Full_Output_fileName = $Full_Output_fileName"

       echo "Starting ease_grid_to_obs for date $Date_output_file"

       # Convert this days worth of obs to a single file 'obs_seq.out'
       # must rename this file to have the date expected by CLM/DART
       ./ease_grid_to_obs

       # link the new observation sequence file to the filename
       # expected by "ease_grid_to_obs"
       \mv -v  obs_seq.out ${Full_Output_fileName}

       # move to next loop for the following day
       echo "Finished ease_grid_to_obs for date $Date_output_file"
   else
      echo "WARNING : No observation files for $YEAR $DOY ... AKA ..."
      echo "WARNING : ${RAW_OBS_DIR}/${YEAR}/${FileBaseDate}*[HV]"
   endif

   @ DOY ++
   @ IDAY ++

end

echo "Finish all data processing"

exit 0

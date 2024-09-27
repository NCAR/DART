#!/bin/ksh -aeux
#
   export PAST_YYYY=2014
   export PAST_MM=07
   export PAST_DD=14
   export YYYY=2014
   export MM=07
   export DD=15
#
# SET OMI PARAMETERS
   export OMI_FILE_PRE=OMI-Aura_L2-OMTO3_
   export OMI_FILE_EXT=.he5
   export OMI_OUTFILE=\'OMI_O3_2014071400.dat\'
#
# SET OBS_WINDOW
   export BIN_BEG=21
   export BIN_END=3
#
# SET OMI INPUT DATA DIR
   export FLG=0
   if [[ ${BIN_END} -eq 3 ]]; then
      export FLG=1
      export BIN_END=24
   fi
   export INFILE=\'/projects/mizzi/TEMPO_FILES/OMI_DATA/OMI_TOO3_20140714_20140717/${OMI_FILE_PRE}${PAST_YYYY}m${PAST_MM}${PAST_DD}t\'
   export OUTFILE=\'TEMP_FILE.dat\'
#
. /etc/profile.d/lmod.sh
module load matlab
matlab -nosplash -nodesktop -r 'omi_o3_extract_test(${INFILE},${OUTFILE},${YYYY},${MM},${DD},${BIN_BEG},${BIN_END},0.,360.,-90.,90.)'
exit


#!/bin/ksh -aeux
export PROJ_NUMBER_ACD=P19010000
export MODIS_INDIR="'"/glade/p/acd/mizzi/AVE_TEST_DATA/MODIS/2008060100"'"
export MODIS_OUTFILE="'"/glade/scratch/mizzi/MODIS/modis_aod_ascii_2008060106_idl"'"
export START_YEAR=2008
export START_MONTH=6
export START_DAY=1
export BIN_BEG=3
export BIN_END=9
export NNL_MIN_LAT=27.
export NNL_MAX_LAT=48.
export NNL_MIN_LON=-132.
export NNL_MAX_LON=-94.
#
# COPY EXECUTABLE
   rm -rf job.ksh
   rm -rf idl_*.err
   rm -rf idl_*.out
   touch job.ksh
   RANDOM=$$
   export JOBRND=idl_${RANDOM}
   cat <<EOFF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER_ACD}
#BSUB -n 1
#BSUB -J ${JOBRND}
#BSUB -o ${JOBRND}.out
#BSUB -e ${JOBRND}.err
#BSUB -W 00:3
#BSUB -q geyser
#
idl << EOF
.compile modis_extract_hdf.pro
modis_extract_hdf, ${MODIS_INDIR}, ${MODIS_OUTFILE}, ${START_YEAR}, ${START_MONTH}, ${START_DAY}, ${BIN_BEG}, ${BIN_END}, ${NNL_MIN_LON}, ${NNL_MAX_LON}, ${NNL_MIN_LAT}, ${NNL_MAX_LAT}
exit
EOF
EOFF
   bsub -K < job.ksh
#

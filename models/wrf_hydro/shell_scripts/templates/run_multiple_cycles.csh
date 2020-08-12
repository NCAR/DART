#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#==========================================================================
# PBS directives                        qsub test_batch.csh
#                                       qstat -u $USER
#                                       qdel <jobnumber>
#PBS -N 6mile_filter
#PBS -e 6mile_filter.stderr
#PBS -o 6mile_filter.stdout
#PBS -l walltime=00:02:00
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -A NRAL0017
#PBS -q premium
#PBS -m abe
#PBS -M thoar@ucar.edu

#==========================================================================

set   MOVE = '/bin/mv -v'
set   COPY = '/bin/cp -vp'
set   LINK = '/bin/ln -vs'
set REMOVE = '/bin/rm -rf'

#==========================================================================
#
# the EXPERIMENT DIRECTORY, a job name, and how many jobs to submit.

set NCYCLES = 10
set RUNDIR = EXPERIMENT_DIRECTORY
set RUNDIR = /glade/scratch/thoar/wrf_hydro/DA_w_params_i2

set YEAR = 2013
set MONTH = 6
set DAY = 27
set HOUR = 22
set MINUTE = 0

set DATESTRING = `printf %04d-%02d-%02d_%02d:%02d $YEAR $MONTH $DAY $HOUR $MINUTE`
set HRLDATE = `printf %04d%02d%02d%02d $YEAR $MONTH $DAY $HOUR`

cd $RUNDIR || exit 1

#==========================================================================
# STEP 1: OBSERVATION
# get DATE_TIME so we can link to correct precomputed identity obs file 
#==========================================================================

${REMOVE} obs_seq.daily

set OBSDAY = `printf %04d%02d%02d $YEAR $MONTH $DAY`
set OBSERVATION_DIR = /glade/work/gharamti/wrfhydro_dart/sixmile/obs_seqs/DA_w_params_i2/USGS
set OBSFILE = ${OBSERVATION_DIR}/obs_seq.${OBSDAY}

# obs_sequence_tool read    an input  file 'obs_seq.daily' 
# obs_sequence_tool creates an output file 'obs_seq.window' 

${COPY} input.nml input.nml.backup

if ( -e    ${OBSFILE} ) then
   ${LINK} ${OBSFILE} obs_seq.daily

   # Chop obs_seq.daily down to desired hour  -30 min+1 sec, +30 min

   set timestring = `printf %04d%02d%02d%02d $YEAR $MONTH $DAY $HOUR`
   set start_time = `echo ${timestring} -1799s -g | ./advance_time` 
   set  last_time = `echo ${timestring} +1800s -g | ./advance_time`

   sed -e "s/first_obs_days.*/first_obs_days = $start_time[1]/"  \
       -e "s/first_obs_seconds.*/first_obs_seconds = $start_time[2]/" \
       -e "s/last_obs_days.*/last_obs_days = $last_time[1]/" \
       -e "s/last_obs_seconds.*/last_obs_seconds = $last_time[2]/" \
       -i input.nml || exit 2

   ./obs_sequence_tool
    
else
   echo "ERROR: No observation file '"${OBSFILE}"' ... exiting."
   exit 3
endif

#==========================================================================
# STEP 2: INFLATION
# IF we are doing inflation, we must take the output inflation files from
# the previous cycle and rename them for input to the current cycle.
#==========================================================================

# We have to potentially deal with files like:
# output_priorinf_mean_d01.${OLDDATESTRING}.nc
# output_priorinf_mean_d02.${OLDDATESTRING}.nc
# output_priorinf_mean_d03.${OLDDATESTRING}.nc
# output_priorinf_sd_d01.${OLDDATESTRING}.nc
# output_priorinf_sd_d02.${OLDDATESTRING}.nc
# output_priorinf_sd_d03.${OLDDATESTRING}.nc
# I am not going to worry about posterior inflation files.

# Should the setup script just create input inflation files so we don't 
# have to screw with changing the namelist after the first execution
# (which traditionally reads from the namelist, not the file)

# If the file exists, just link to the new expected name.
# the expected names have a _d0? inserted before the file extension
# if there are multiple domains.
# If the file does not exist, filter will die and issue a very explicit
# death message.

${REMOVE} input_priorinf_mean*.nc input_priorinf_sd*.nc

set n_domain = `grep domain_order input.nml | cut -d'=' -f2 | tr -d ' ' | tr '"' "'" | tr "," '\n' | egrep -v '^$' | wc -l`

if ( $n_domain == 1 ) then
    set the_domains = 'single_domain'
else
    set the_domains = "_d01 _d02 _d03"
endif

foreach DOMAIN ( `echo $the_domains`  )

   if ( $DOMAIN == "single_domain" ) then
       set DOMAIN = ''
   endif

   # Checking for a prior inflation mean file from the previous assimilation.

   (ls -rt1 output_priorinf_mean${DOMAIN}.* | tail -n 1 >! latestfile) > & /dev/null
   set nfiles = `cat latestfile | wc -l`

   if ( $nfiles > 0 ) then
      set latest = `cat latestfile`
      ${LINK} $latest input_priorinf_mean${DOMAIN}.nc
   endif

   # Checking for a prior inflation sd file from the previous assimilation.

   (ls -rt1 output_priorinf_sd${DOMAIN}.* | tail -n 1 >! latestfile) > & /dev/null
   set nfiles = `cat latestfile | wc -l`

   if ( $nfiles > 0 ) then
      set latest = `cat latestfile`
      ${LINK} $latest input_priorinf_sd${DOMAIN}.nc
   endif

end

${REMOVE} latestfile

#==========================================================================
# STEP 3: Assimilate.
# Run DART on the ensemble of new states.
# Collect all the RESTARTs for each domain into a list of input files. 
# The io module will error out if the input file list is too short 
# which helps make sure all instances advanced successfully.
# Our strategy is that DART (filter) will modify these files in-place.
# If you need to save a copy, do so now, or set one of the DART 
# 'stages_to_write' to 'input' and 'num_output_state_members = ens_size'
# and 'output_members = .true.'. This will write _minimal_ netCDF files
# with whatever is in the DART state. You could take these variables and 
# insert them into a 'full' restart file and run ...
#==========================================================================

# Clean up from any previous execution
${REMOVE} dart_log.out dart_log.nml
${REMOVE} lsm_file_list.txt hydro_file_list.txt  param_file_list.txt

@ ens_size = 0
foreach MEMBER ( member_* )
   ls -rt1 $MEMBER/RESTART.*DOMAIN*   | tail -n 1 >> lsm_file_list.txt
   ls -rt1 $MEMBER/HYDRO_RST.*DOMAIN* | tail -n 1 >> hydro_file_list.txt
   ls -rt1 $MEMBER/param*.nc          | tail -n 1 >> param_file_list.txt
   @ ens_size ++
end

# If there are no files for that domain ... just remove the (empty) file.

if (`cat lsm_file_list.txt | wc -l` != $ens_size) then
   ${REMOVE} lsm_file_list.txt
endif

if (`cat hydro_file_list.txt | wc -l` != $ens_size) then
   ${REMOVE} hydro_file_list.txt
endif

if (`cat param_file_list.txt | wc -l` != $ens_size) then
   ${REMOVE} param_file_list.txt
endif

# Perform the assimilation.

mpiexec_mpt ./filter || exit 4

# Tag the output with the valid time of the model state.
# TODO could move each ensemble-member file to the respective member dir.

foreach FILE ( input_*mean.nc      input_*sd.nc \
            forecast_*mean.nc   forecast_*sd.nc  forecast_member_????.nc \
            preassim_*mean.nc   preassim_*sd.nc  preassim_member_????.nc \
           postassim_*mean.nc  postassim_*sd.nc postassim_member_????.nc \
            analysis_*mean.nc   analysis_*sd.nc  analysis_member_????.nc \
              output_*mean.nc     output_*sd.nc \
              output_*mean_d0?.nc output_*sd_d0?.nc )

   if (  -e $FILE ) then
      set FEXT  = $FILE:e
      set FBASE = $FILE:r
      ${MOVE} $FILE ${FBASE}.${DATESTRING}.${FEXT}
   else
      echo "$FILE does not exist, no need to take action."
   endif
end

# Tag the DART observation file with the valid time of the model state.

${MOVE} obs_seq.final obs_seq.final.${DATESTRING}

#==========================================================================
# STEP 4: If necessary, submit the job to advance the ensemble.
# look at gps/shell_scripts/multi_parallel.batch
#==========================================================================

set imem = 0

while ( $imem < $ens_size )

   @ mynum = $imem + 1
   echo "advancing member $mynum of $ens_size"
   set MYDIR = `printf member_%03d $imem`
   cd $RUNDIR/$MYDIR

   set hydroFILE = "./HYDRO_RST.${DATESTRING}_DOMAIN1"
   set hrldasFILE = "./RESTART.${HRLDATE}_DOMAIN1"

   echo "set hydroFILE = ${hydroFILE}"
   echo "set hrldasFILE = ${hrldasFILE}"

   sed -e "s#restart_file.*#restart_file = '${hydroFILE}'#" \
          -i hydro.namelist || exit 5

   sed -e "s/START_YEAR.*/START_YEAR = ${YEAR}/" \
       -e "s/START_MONTH.*/START_MONTH = ${MONTH}/" \
       -e "s/START_DAY.*/START_DAY = ${DAY}/" \
       -e "s/START_HOUR.*/START_HOUR = ${HOUR}/" \
       -e "s/START_MIN.*/START_MIN = ${MINUTE}/" \
       -e "s#RESTART_FILENAME_REQUESTED.*#RESTART_FILENAME_REQUESTED = '${hrldasFILE}'#" \
       -i namelist.hrldas || exit 6
 
   mpiexec_mpt ./wrf_hydro.exe || exit 7

   @ imem++

end

exit 0


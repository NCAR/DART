#!/bin/bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms.

#   driver.sh - script that performs assimilation cycling
#
#   run script as: ./driver.sh 2024051906 param.sh >& run.out &# Run the assimilation experiment

main() {

if [[ $# < 2 ]]; then 
   echo "Please enter a date {yyyymmddhh} and the paramfile"
   exit 0
fi

datea=$1                     # starting date
paramfile=$(readlink -f $2)  # param file

# Target date: last assimilation time
datefnl=2024051912

# Max WRF re-runs
MAX_WRF_RERUNS=1

# 1. Make sure the input arguments are valid
validate_input

# 2. Main assimilation loop
touch $RUN_DIR/cycle_started_$datea

while true; do

   # Previous assimilation time (usually the date before 6 hours)
   datep=$(echo $datea -$ASSIM_INT_HOURS | $DART_DIR/models/wrf/work/advance_time)

   configure_da_init
   prepare_files_for_da

   run_filter
   post_da_diags

   integrate_wrf 
   cycle_cleanup 

   # Are we done?
   if [[ $datea == $datefnl ]]; then
      echo Reached the final date 
      echo Script exiting normally
      exit 0
   fi

   # Next assimilation time
   datea=$(echo $datea $ASSIM_INT_HOURS | $DART_DIR/models/wrf/work/advance_time)

done # Main Loop


}



###############
## FUNCTIONS ##
###############

# ------------------------------
# Manage input and namelist options

validate_input() {

source $paramfile; cd $RUN_DIR

$COPY $TEMPLATE_DIR/namelist.input.meso  namelist.wps.check

hybrid_opt=$(grep  "hybrid_opt"  namelist.wps.check | cut -d"=" -f2 | xargs)
use_theta_m=$(grep "use_theta_m" namelist.wps.check | cut -d"=" -f2 | xargs)

$REMOVE namelist.wps.check

if [[ $hybrid_opt != 0 ]]; then
   echo "ERROR: WRF-DART must use terrain following coordinates, hybrid_opt must be 0"
   exit 1
fi
if [[ $use_theta_m != 0 ]]; then
   echo "ERROR: WRF-DART must assign WRF variable THM as T, use_theta_m must be 0"
   exit 1
fi

}


# ------------------------------
# Manage time and output directories

configure_da_init() {

# Clean any files from previous runs
$REMOVE $RUN_DIR/cycle_finished_* \
        $RUN_DIR/prev_cycle_done  || true

# Check if the output directory exists
if [[ ! -d $OUTPUT_DIR/$datea ]]; then
   $REMOVE $RUN_DIR/ABORT_RETRO
   echo "Exiting because output directory does not exist."
   exit 2
fi   

read -r -a gdate  < <(echo $datea                0 -g | $DART_DIR/models/wrf/work/advance_time)
read -r -a gdatef < <(echo $datea $ASSIM_INT_HOURS -g | $DART_DIR/models/wrf/work/advance_time)

wdate=$(echo $datea 0 -w | $DART_DIR/models/wrf/work/advance_time)
hh=${datea:8:2}

echo "Ready to check inputs"

for (( dn=1; dn<=NUM_DOMAINS; dn++ )); do
    
    dchar=$(echo $dn + 100 | bc | cut -b2-3)
    
    # Check if the boundary/initial WRF files 
    # and observations are available
    for infile in \
        wrfinput_d${dchar}_${gdate[0]}_${gdate[1]}_mean   \
        wrfinput_d${dchar}_${gdatef[0]}_${gdatef[1]}_mean \
        wrfbdy_d01_${gdatef[0]}_${gdatef[1]}_mean         \
        obs_seq.out; do

        if [[ ! -e $OUTPUT_DIR/$datea/$infile ]]; then
           echo  "$OUTPUT_DIR/$datea/$infile is missing! Stopping the system"
           touch ABORT_RETRO
           exit 3
        fi
    done

done

}


# ------------------------------
# Extract variables for DA

prepare_files_for_da() {

# Check if the prior ensemble is available 
# and do some cleanup in the rundir
for (( dn=1; dn<=NUM_DOMAINS; dn++ )); do
  
    dchar=$(echo $dn + 100 | bc | cut -b2-3)

    for (( n=1; n<=NUM_ENS; n++ )); do   
           
        ensstring=$(echo $n + 10000 | bc | cut -b2-5)
   
        if [[ -e "$OUTPUT_DIR/$datep/PRIORS/prior_d$dchar.$ensstring" ]]; then
           if (( dn == 1 )) && [[ -d $RUN_DIR/advance_temp$n ]]; then
               $REMOVE $RUN_DIR/advance_temp$n
           fi
           mkdir -p $RUN_DIR/advance_temp$n
           $LINK    $OUTPUT_DIR/$datea/wrfinput_d${dchar}_${gdate[0]}_${gdate[1]}_mean \
                    $RUN_DIR/advance_temp$n/wrfinput_d$dchar
        else
           echo "$OUTPUT_DIR/$datep/PRIORS/prior_d$dchar.$ensstring is missing! Stopping the system"
           touch ABORT_RETRO
           exit 4
        fi
    
    done

done

# Extract WRF variables to update (for each member) 
# This will call another script "prep_ic.sh"
for (( n=1; n<=NUM_ENS; n++ )); do
    $SHELL_SCRIPTS_DIR/prep_ic.sh $n $datep $NUM_DOMAINS $paramfile
done

# Create WRF template files
for (( dn=1; dn<=NUM_DOMAINS; dn++ )); do
    dchar=$(echo $dn + 100 | bc | cut -b2-3)
    $COPY $OUTPUT_DIR/$datea/wrfinput_d${dchar}_${gdate[0]}_${gdate[1]}_mean wrfinput_d$dchar
done

# Move inflation files 
if [[ $ADAPTIVE_INFLATION == 1 ]]; then

   # Create the home for inflation and future state space diagnostic files
   mkdir -p $RUN_DIR/Inflation_input $RUN_DIR/Output

   for inf_file in $OUTPUT_DIR/$datep/Inflation_input/*.nc ; do
       $LINK $inf_file $RUN_DIR 
   done  
fi     

$LINK   $OUTPUT_DIR/$datea/obs_seq.out $RUN_DIR
$REMOVE $RUN_DIR/WRF
$REMOVE $RUN_DIR/prev_cycle_done
$LINK   $OUTPUT_DIR/$datea $RUN_DIR/WRF

mkdir -p $OUTPUT_DIR/$datea/logs

}


# ------------------------------
# Perform the DA to generate the analysis

run_filter() {

$REMOVE script.sed || true

assimilate_job=$RUN_DIR/assimilate.sh

cat > $assimilate_job <<EOF
    #!/bin/bash
    #=================================================================
    #PBS -N assimilate_$datea
    #PBS -j oe
    #PBS -A $COMPUTER_CHARGE_ACCOUNT
    #PBS -l walltime=$FILTER_TIME
    #PBS -q $FILTER_QUEUE
    #PBS -l job_priority=$FILTER_PRIORITY
    #PBS -m ae
    #PBS -M $EMAIL
    #PBS -k eod
    #PBS -l select=$FILTER_NODES:ncpus=$FILTER_PROCS:mpiprocs=$FILTER_MPI
    #=================================================================

    $SHELL_SCRIPTS_DIR/assimilate.sh $datea $paramfile
EOF

chmod +x $assimilate_job
qsub $assimilate_job

filter_thresh_min=$(echo $FILTER_TIME | cut -b3-4)
filter_thresh_hr=$(echo $FILTER_TIME | cut -b1-1)
filter_thresh=$(( 10#${filter_thresh_min} * 60 + 10#${filter_thresh_hr} * 3600 ))

while [[ ! -e filter_done ]]; do

   # Check the timing.  If it took longer than the time allocated, abort.
   if [[ -e filter_started ]]; then

      start_time=$(head -1 filter_started)
      end_time=$(date +%s)
      total_time=$(( end_time - start_time ))

      if (( total_time > filter_thresh )); then
         # If the job needs to be aborted ... we need to qdel the hanging job
         echo "Time exceeded the maximum allowable time. Exiting."
         touch ABORT_RETRO
         $REMOVE filter_started
         exit 5
      fi

   fi
   sleep 10

done

}


# ------------------------------
# Clean up and move files

post_da_diags() {

echo "Filter is done, cleaning up"

$REMOVE $RUN_DIR/filter_started      \
        $RUN_DIR/filter_done         \
        $RUN_DIR/obs_seq.out         \
        $RUN_DIR/ic_*_ready          \
        $RUN_DIR/postassim_priorinf* \
        $RUN_DIR/preassim_priorinf*  || true

if [[ -e assimilate.sh ]]; then
   $REMOVE $RUN_DIR/assimilate.sh
fi

mkdir -p ${OUTPUT_DIR}/${datea}/Inflation_input $OUTPUT_DIR/$datea/WRFIN $OUTPUT_DIR/$datea/PRIORS          

# Variables to extract for the increment
extract_str=""
if declare -p increment_vars_a &>/dev/null; then
   for v in ${increment_vars_a[@]}; do
       if [[ -z ${extract_str} ]]; then
          extract_str=${v}
       else
          extract_str="${extract_str},${v}"
       fi
   done
fi

# Increments
for (( dn=1; dn<=NUM_DOMAINS; dn++ )); do
    
    dchar=$(echo $dn + 100 | bc | cut -b2-3)

    # Could also read in the requested stages from the input.nml
    if (( NUM_DOMAINS == 1 )); then
       ncdiff -F -O -v    $extract_str postassim_mean.nc preassim_mean.nc analysis_increment.nc
       ncks   -F -O -x -v $extract_str postassim_mean.nc static_data.nc
       ncks   -A                       static_data.nc analysis_increment.nc
    else
       ncdiff -F -O -v    $extract_str postassim_mean_d$dchar.nc preassim_mean_d$dchar.nc analysis_increment_d$dchar.nc
       ncks   -F -O -x -v $extract_str postassim_mean_d$dchar.nc static_data_d$dchar.nc
       ncks   -A                       static_data_d$dchar.nc analysis_increment_d$dchar.nc
    fi

done 

# Move diagnostic and obs_seq.final data to storage directories
files=( postassim_*.nc preassim_*.nc
        analysis_*.nc  forecast_*.nc
        output_m*.nc   output_s*.nc 
	obs_seq.final
        	
      )
    
for FILE in ${files[@]}; do    
    if [[ -s $FILE ]]; then
       $MOVE $FILE $OUTPUT_DIR/$datea
    fi
done

# Move inflation files
if [[ $ADAPTIVE_INFLATION == 1 ]]; then

   for (( dn=1; dn<=NUM_DOMAINS; dn++ )); do
       
       dchar=$(echo $dn + 100 | bc | cut -b2-3)
       
       if (( NUM_DOMAINS ==  1 )); then
	  if [[ -e output_priorinf_mean.nc ]]; then
             echo "Detected output_priorinf_mean.nc during post_da_diags moves"   
          else
             echo "Error: Missing output_priorinf_mean.nc during post_da_diags moves"
             exit	     
          fi     
          old_file=( input_postinf_mean.nc  input_postinf_sd.nc  input_priorinf_mean.nc  input_priorinf_sd.nc )
          new_file=( output_postinf_mean.nc output_postinf_sd.nc output_priorinf_mean.nc output_priorinf_sd.nc ) 
       else
          if [[ -e output_priorinf_mean_d${dchar}.nc ]]; then
	     echo "Detected output_priorinf_mean_d${dchar}.nc  during post_da_diags moves"   
	  else
             echo "Error: Missing output_priorinf_mean_d${dchar}.nc  during post_da_diags moves"
             exit
          fi
	  old_file=( input_postinf_mean_d$dchar.nc   input_postinf_sd_d$dchar.nc
                     input_priorinf_mean_d$dchar.nc  input_priorinf_sd_d$dchar.nc  )
          new_file=( output_postinf_mean_d$dchar.nc  output_postinf_sd_d$dchar.nc
                     output_priorinf_mean_d$dchar.nc output_priorinf_sd_d$dchar.nc )
       fi
       nfiles=${#new_file[@]}
       for (( i=0; i<nfiles; i++ )); do        
           if [[ -e ${new_file[$i]} && -s ${new_file[$i]} ]]; then
              $MOVE ${new_file[$i]} $OUTPUT_DIR/$datea/Inflation_input/${old_file[$i]}
           fi
       done

   done  # loop through domains

   echo "Past the inflation file moves."
fi

# Obs-space diagnostics
if [[ -e obs_diag.log ]]; then
   $REMOVE obs_diag.log
fi

$SHELL_SCRIPTS_DIR/diagnostics_obs.sh $datea $paramfile >& ${RUN_DIR}/obs_diag.log  

}



#-------------------------------
# Run the ensemble forward in time 

integrate_wrf() {

# Clean any previous log files
if [[ -e $RUN_DIR/start_member_1 ]]; then
   $REMOVE $RUN_DIR/start_member_* $RUN_DIR/done_member_* 
fi

# Submission scripts
for (( n=1; n<=NUM_ENS; n++ )); do

    jobfile=assim_advance_mem$n.sh
   
    cat > $jobfile <<EOF
        #!/bin/bash
        #=================================================================
        #PBS -N assim_advance_$n
        #PBS -j oe
        #PBS -A $COMPUTER_CHARGE_ACCOUNT
        #PBS -l walltime=$ADVANCE_TIME
        #PBS -q $ADVANCE_QUEUE
        #PBS -l job_priority=$ADVANCE_PRIORITY
        #PBS -m a
        #PBS -M $EMAIL
        #PBS -k eod
        #PBS -l select=$ADVANCE_NODES:ncpus=$ADVANCE_PROCS:mpiprocs=$ADVANCE_MPI
        #=================================================================

        set -uo pipefail

        cd $RUN_DIR

        # Run the actual advance script (no scheduler headers inside it)
        $SHELL_SCRIPTS_DIR/assim_advance.sh $datea $n $paramfile
EOF
   
    chmod +x $jobfile
    qsub $jobfile

done

# Make sure all members have advanced
advance_thresh_min=$(echo $ADVANCE_TIME | cut -b3-4)
advance_thresh_hr=$(echo $ADVANCE_TIME  | cut -b1-1)
advance_thresh=$(( 10#${advance_thresh_min} * 60 + 10#${advance_thresh_hr} * 3600 ))

for (( n=1; n<=NUM_ENS; n++ )); do

    ensstring=$(echo $n + 10000 | bc | cut -b2-5)
    keep_trying=true
    max_retry=1

    while [[ $keep_trying == true ]]; do

        # Wait for the script to start
        while [[ ! -e $RUN_DIR/start_member_$n ]]; do
  
           if [[ $(qstat -wa | grep -c assim_advance_$n || true) -eq 0 ]]; then
              echo "Warning, detected that assim_advance_$n is missing from the queue"
              echo "If this warning leads to  missing output from ensemble $n"
              echo "consider enabling the qsub command within keep_trying while statement in driver.sh"
           fi
           sleep 5

        done
       
        start_time=$(head -1 start_member_$n)
        echo "Member $n has started.  Start time $start_time"

        # Wait for the output file
        while true; do

           current_time="$(date +%s)"
           length_time=$(( current_time - start_time ))

           if [[ -e $RUN_DIR/done_member_$n ]]; then

              #  If the output file already exists, move on
              keep_trying=false
              break
       
           elif (( length_time > advance_thresh )); then

              #  If WRF member has failed 1 resubmission attempt, immediately stop driver.sh
              if (( max_retry > MAX_WRF_RERUNS )); then

                 echo "Stopping the driver.sh script! The WRF ensemble member $n"
                 echo "has exceeded the maximum resubmission attempts (1) without completing."
                 echo "This typically means the WRF integration has failed."
                 echo "Check your BASE_DIR/rundir/advance_temp$n directory and locate"
                 echo "the WRF rsl.out.0000 or rsl.error.0000 log files for further information."
                 echo "If applicable, check the DART analysis_increment.nc from previous assimilation step"
                 exit 6

              fi

              # The WRF job did not complete. Resubmit to queue
              $REMOVE "start_member_$n"
              echo "Did not find the member done file, WRF run did not complete"
              echo "Attempting resubmission $max_retry"
              (( max_retry++ ))

              qsub "assim_advance_mem$n.sh"
              sleep 5 
           fi
           sleep 15 
           
        done 
    done

    for (( dn=1; dn<=NUM_DOMAINS; dn++ )); do

        dchar=$(echo $dn + 100 | bc | cut -b2-3)
        
        echo "Moving $n $ensstring for domain $dn"

        # Move log files and clean up
        if (( dn == 1 )); then
           $MOVE   $RUN_DIR/assim_advance_$n.o*                $OUTPUT_DIR/$datea/logs
           $MOVE   WRFOUT/wrf.out_${gdatef[0]}_${gdatef[1]}_$n $OUTPUT_DIR/$datea/logs
           $REMOVE start_member_$n done_member_$n

           if [[ -e assim_advance_mem$n.sh ]]; then
              $REMOVE assim_advance_mem$n.sh
           fi

           pert=$(cat $RUN_DIR/advance_temp$n/mem${n}_pert_bank_num)
           echo "Member $n uses perturbation bank ensemble member $pert" >> "$OUTPUT_DIR/$datea/pert_bank_members.txt"
        fi

        $MOVE   WRFIN/wrfinput_d${dchar}_$n.gz    $OUTPUT_DIR/$datea/WRFIN
        $MOVE   $RUN_DIR/prior_d$dchar.$ensstring $OUTPUT_DIR/$datea/PRIORS
        $REMOVE filter_restart_d$dchar.$ensstring
      
    done 
done

}



#-------------------------------
# Cleanup and finish the cycle 

cycle_cleanup() {

if [[ -e obs_prep.log ]]; then
   $REMOVE obs_prep.log
fi

# Move data to storage directory
$COPY input.nml $OUTPUT_DIR/$datea
$MOVE $RUN_DIR/assimilate*.o* \
      $RUN_DIR/dart_log.*     \
      $RUN_DIR/*.log          \
      $RUN_DIR/assim*.o*      \
      $OUTPUT_DIR/$datea/logs || true

# Remove temporary files from both the run directory and old storage directories
$REMOVE $OUTPUT_DIR/$datep/wrfinput_d*_mean $RUN_DIR/wrfinput_d* $RUN_DIR/WRF

# Prep data for archive
cd $OUTPUT_DIR/$datea

gzip  -f wrfinput_d*_${gdate[0]}_${gdate[1]}_mean wrfinput_d*_${gdatef[0]}_${gdatef[1]}_mean wrfbdy_d*_mean
tar -cvf retro.tar wrfin*.gz wrfbdy_d*.gz

if [[ $ADAPTIVE_INFLATION == 1 ]]; then
   tar -rvf dart_data.tar obs_seq.out obs_seq.final Inflation_input logs input.nml
else
   tar -rvf dart_data.tar obs_seq.out obs_seq.final logs input.nml
fi

$REMOVE wrfinput_d*_${gdate[0]}_${gdate[1]}_mean.gz wrfbdy_d*.gz
gunzip -f wrfinput_d*_${gdatef[0]}_${gdatef[1]}_mean.gz

cd $RUN_DIR

$REMOVE $RUN_DIR/input_priorinf_* || true 
$REMOVE $RUN_DIR/static_data*     || true

touch prev_cycle_done
touch $RUN_DIR/cycle_finished_$datea

if [[ -e cycle_started_$datea ]]; then
   rm -f $RUN_DIR/cycle_started_$datea
fi

}

#-------------------------------

main "$@"

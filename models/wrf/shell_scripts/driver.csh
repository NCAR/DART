#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#   driver.csh - script that is the driver for the
#                            CONUS analysis system
#                            MODIFIED for new DART direct
#                            file access
#
#      provide an input argument of the first
#      analysis time in yyyymmddhh format.
#
#   Created May 2009, Ryan Torn, U. Albany
#   Modified by G. Romine to run realtime cases 2011-18
#
########################################################################
#   run as: nohup csh driver.csh 2017042706 param.csh >& run.log &
########################################################################
# Set the correct values here
set paramfile = `readlink -f ${2}` # Get absolute path for param.csh from command line arg
set datefnl   =  2017042712 # target date   YYYYMMDDHH  # set this appropriately #%%%#
########################################################################
# Likely do not need to change anything below
########################################################################

source $paramfile

echo `uname -a`
cd ${RUN_DIR}

#  First determine the appropriate analysis date

if ( $#argv > 0 ) then
  set datea   = ${1} # starting date
  setenv restore 1   # set the restore variable
  echo 'starting a restore'
else
  echo "please enter a date: yyyymmddhh"
  exit
endif

touch $RUN_DIR/cycle_started_${datea}

while ( 1 == 1 )

   if ( ! -d ${OUTPUT_DIR}/${datea} && $restore == 1 ) then        	
      ${REMOVE} ${RUN_DIR}/ABORT_RETRO
      echo 'exiting because output directory does not exist and this is a restore'
      exit
   endif

   set datep  = `echo $datea -${ASSIM_INT_HOURS}   | ${DART_DIR}/models/wrf/work/advance_time`
   set gdate  = `echo $datea 0 -g                  | ${DART_DIR}/models/wrf/work/advance_time`
   set gdatef = `echo $datea ${ASSIM_INT_HOURS} -g | ${DART_DIR}/models/wrf/work/advance_time`
   set wdate  = `echo $datea 0 -w                  | ${DART_DIR}/models/wrf/work/advance_time`
   set hh     = `echo $datea | cut -b9-10`

   echo 'ready to check inputs'
   set domains = $NUM_DOMAINS   # from the param file
   #  Check to make sure all input data exists
   if  ( $domains == 1 ) then
      foreach infile ( wrfinput_d01_${gdate[1]}_${gdate[2]}_mean \
                      wrfinput_d01_${gdatef[1]}_${gdatef[2]}_mean \
                        wrfbdy_d01_${gdatef[1]}_${gdatef[2]}_mean obs_seq.out )

         if ( ! -e ${OUTPUT_DIR}/${datea}/${infile} ) then
            echo  "${OUTPUT_DIR}/${datea}/${infile} is missing!  Stopping the system"
            touch ABORT_RETRO
            exit 2
         endif
      end
   endif

   #  Clear the advance_temp directory, write in new template file, and
   # overwrite variables with the compact prior netcdf files
   #
   #  NOTE that multiple domains might be present, but only looking for domain 1

   if ( $SUPER_PLATFORM == 'LSF queuing system' ) then
      set ic_queue = caldera
      set logfile = "${RUN_DIR}/ic_gen.log"
      set sub_command = "bsub -q ${ic_queue} -W 00:05 -o ${logfile} -n 1 -P ${COMPUTER_CHARGE_ACCOUNT}"
   else if ( $SUPER_PLATFORM == 'derecho' ) then
      set ic_queue = "main"
      set sub_command = "qsub -l select=1:ncpus=128:mpiprocs=128:mem=5GB -l walltime=00:03:00 -q ${ic_queue} -A ${COMPUTER_CHARGE_ACCOUNT} -j oe -k eod -N icgen "
   endif

   echo "this platform is $SUPER_PLATFORM and the job submission command is $sub_command"

   set dn = 1
   while ( $dn <= $domains )
      set dchar = `echo $dn + 100 | bc | cut -b2-3`
      set n = 1
      while ( $n <= $NUM_ENS )
         set ensstring = `echo $n + 10000 | bc | cut -b2-5`
         if ( -e ${OUTPUT_DIR}/${datep}/PRIORS/prior_d${dchar}.${ensstring} ) then

            if ( $dn == 1 &&  -d ${RUN_DIR}/advance_temp${n} )  ${REMOVE} ${RUN_DIR}/advance_temp${n}

            mkdir -p ${RUN_DIR}/advance_temp${n}
            ${LINK} ${OUTPUT_DIR}/${datea}/wrfinput_d${dchar}_${gdate[1]}_${gdate[2]}_mean \
                       ${RUN_DIR}/advance_temp${n}/wrfinput_d${dchar}
         else
            echo "${OUTPUT_DIR}/${datep}/PRIORS/prior_d${dchar}.${ensstring} is missing! Stopping the system"
            touch ABORT_RETRO
            exit 3
         endif
         @ n++
      end  # loop through ensemble members
      @ dn++
   end   # loop through domains

   # Fire off a bunch of small jobs to create the initial conditions for the short model forecast.
   # the prep_ic.csh script creates a file "${RUN_DIR}/ic_d${dchar}_${n}_ready" to signal a
   # successful completion.
   # NOTE : Submit commands here are system specific and work for this tutorial, users may want/need to change
   #        for their system and/or production.

   set n = 1
   while ( $n <= $NUM_ENS )
      if ( $SUPER_PLATFORM == 'derecho' ) then   # can't pass along arguments in the same way
         $sub_command -v mem_num=${n},date=${datep},domain=${domains},paramf=${paramfile} ${SHELL_SCRIPTS_DIR}/prep_ic.csh
      else
         $sub_command " ${SHELL_SCRIPTS_DIR}/prep_ic.csh ${n} ${datep} ${dn} ${paramfile} "
      endif
      @ n++
   end  # loop through ensemble members

   # If any of the queued jobs has not completed in 5 minutes, run them manually
   # cleanup any failed stuffs
   # NOTE : No automated cleanup for queued jobs. User may want to add system specific monitoring.
   set dn = 1
   while ( $dn <= $domains )
      set dchar = `echo $dn + 100 | bc | cut -b2-3`
      set n = 1
      set loop = 1
      while ( $n <= $NUM_ENS )
         if (  -e    ${RUN_DIR}/ic_d${dchar}_${n}_ready) then
            ${REMOVE} ${RUN_DIR}/ic_d${dchar}_${n}_ready
            @ n++
            set loop = 1
         else
            echo "waiting for ic member $n in domain $dn"
            sleep 5
            @ loop++
            if ( $loop > 60 ) then    # wait 5 minutes for the ic file to be ready, else run manually
               echo "gave up on ic member $n - redo"
               ${SHELL_SCRIPTS_DIR}/prep_ic.csh ${n} ${datep} ${dn} ${paramfile}
               # If manual execution of script, shouldn't queued job be killed?
            endif
         endif
      end
      @ dn++
   end   # loop through domains

   mkdir ${OUTPUT_DIR}/${datea}/logs
   ${MOVE}  icgen.o* ${OUTPUT_DIR}/${datea}/logs/

   #  Get wrfinput source information
   ${COPY} ${OUTPUT_DIR}/${datea}/wrfinput_d01_${gdate[1]}_${gdate[2]}_mean wrfinput_d01
   set dn = 1
   while ( $dn <= $domains )

      set dchar = `echo $dn + 100 | bc | cut -b2-3`
      ${COPY} ${OUTPUT_DIR}/${datea}/wrfinput_d${dchar}_${gdate[1]}_${gdate[2]}_mean wrfinput_d${dchar}
      @ dn++

   end

   # Copy the inflation files from the previous time, update for domains
   #TJH ADAPTIVE_INFLATION comes from scripts/param.csh but is disjoint from input.nml

   if ( $ADAPTIVE_INFLATION == 1 ) then
      # Create the home for inflation and future state space diagnostic files
      # Should try to check each file here, but shortcutting for prior (most common) and link them all

      mkdir -p ${RUN_DIR}/{Inflation_input,Output}

      if ( $domains == 1) then
         if ( -e ${OUTPUT_DIR}/${datep}/Inflation_input/input_priorinf_mean.nc ) then

            ${LINK} ${OUTPUT_DIR}/${datep}/Inflation_input/input_priorinf*.nc ${RUN_DIR}/.
            ${LINK} ${OUTPUT_DIR}/${datep}/Inflation_input/input_postinf*.nc ${RUN_DIR}/.

         else

            echo "${OUTPUT_DIR}/${datep}/Inflation_input/input_priorinf_mean.nc file does not exist.  Stopping"
            touch ABORT_RETRO
            exit 3

         endif

      else    # multiple domains so multiple inflation files for each domain
              # TJH this should error out much earlier
         echo "This script doesn't support multiple domains.  Stopping"
         touch ABORT_RETRO
         exit 4

      endif # number of domains check

   endif   # ADAPTIVE_INFLATION file check

   ${LINK} ${OUTPUT_DIR}/${datea}/obs_seq.out .
   ${REMOVE} ${RUN_DIR}/WRF
   ${REMOVE} ${RUN_DIR}/prev_cycle_done
   ${LINK} ${OUTPUT_DIR}/${datea} ${RUN_DIR}/WRF

   #  run filter to generate the analysis
   ${REMOVE} script.sed
   if ( $SUPER_PLATFORM == 'LSF queuing system' ) then

      # This is a most unusual application of 'sed' to insert the batch submission
      # directives into a file. The last backslash '\' before the quote is essential.
      # What happens to the first quote on the next line is beyond me ... TJH.
      # Other places in DART simply have both sets of directives in the template
      # file and sed just replaces singular values.

      echo "2i\"                                                                  >! script.sed
      echo "#==================================================================\" >> script.sed
      echo "#BSUB -J assimilate_${datea}\"                                        >> script.sed
      echo "#BSUB -o assimilate_${datea}.%J.log\"                                 >> script.sed
      echo "#BSUB -P ${COMPUTER_CHARGE_ACCOUNT}\"                                 >> script.sed
      echo "#BSUB -W ${FILTER_TIME}\"                                             >> script.sed
      echo "#BSUB -q ${FILTER_QUEUE}\"                                            >> script.sed
      echo "#BSUB -n ${FILTER_CORES}\"                                            >> script.sed
      echo "#BSUB -x\"                                                            >> script.sed
      echo '#BSUB -R "span[ptile='"${FILTER_PTILE}]"'"\'                          >> script.sed
      echo "#=================================================================="  >> script.sed
      echo 's%${1}%'"${datea}%g"                                                  >> script.sed
      echo 's%${2}%'"${paramfile}%g"                                              >> script.sed
      sed -f script.sed ${SHELL_SCRIPTS_DIR}/assimilate.csh >! assimilate.csh

      if ( $?reservation ) then
         echo "USING RESERVATION," `/contrib/lsf/get_my_rsvid`
         bsub -U `/contrib/lsf/get_my_rsvid` < assimilate.csh
      else
         bsub < assimilate.csh
      endif
      set this_filter_runtime = $FILTER_TIME

   else if ( $SUPER_PLATFORM == 'derecho' ) then

      echo "2i\"                                                                          >! script.sed
      echo "#=================================================================\"          >> script.sed
      echo "#PBS -N assimilate_${datea}\"                                                 >> script.sed
      echo "#PBS -j oe\"                                                                  >> script.sed
      echo "#PBS -A ${COMPUTER_CHARGE_ACCOUNT}\"                                          >> script.sed
      echo "#PBS -l walltime=${FILTER_TIME}\"                                             >> script.sed
      echo "#PBS -q ${FILTER_QUEUE}\"                                                     >> script.sed
      echo "#PBS -l job_priority=${FILTER_PRIORITY}\"                                     >> script.sed
      echo "#PBS -m ae\"                                                                  >> script.sed
      echo "#PBS -M ${EMAIL}\"                                                            >> script.sed
      echo "#PBS -k eod\"                                                                 >> script.sed
      echo "#PBS -l select=${FILTER_NODES}:ncpus=${FILTER_PROCS}:mpiprocs=${FILTER_MPI}\" >> script.sed
      echo "#================================================================="           >> script.sed
      echo 's%${1}%'"${datea}%g"                                                          >> script.sed
      echo 's%${2}%'"${paramfile}%g"                                                      >> script.sed
      sed -f script.sed ${SHELL_SCRIPTS_DIR}/assimilate.csh >! assimilate.csh

      qsub assimilate.csh

      set this_filter_runtime = $FILTER_TIME

   endif

   cd $RUN_DIR   # make sure we are still in the right place

   set filter_thresh = `echo $this_filter_runtime | cut -b3-4`
   @ filter_thresh = `expr $filter_thresh \+ 0` * 60 + `echo $this_filter_runtime | cut -b1-1` * 3600

   while ( ! -e filter_done )

      # Check the timing.  If it took longer than the time allocated, abort.
      if ( -e filter_started ) then

         set start_time = `head -1 filter_started`
         set end_time = `date +%s`

         @ total_time = $end_time - $start_time
         if ( $total_time > $filter_thresh ) then

            # If the job needs to be aborted ... we need to qdel the hanging job

            echo "Time exceeded the maximum allowable time.  Exiting."
            touch ABORT_RETRO
            ${REMOVE} filter_started
            exit 5

         endif

      endif
      sleep 10

   end

   echo "filter is done, cleaning up"

   ${MOVE}  icgen.o* ${OUTPUT_DIR}/${datea}/logs/
   ${REMOVE} ${RUN_DIR}/filter_started  \
             ${RUN_DIR}/filter_done  \
             ${RUN_DIR}/obs_seq.out     \
             ${RUN_DIR}/postassim_priorinf*  \
             ${RUN_DIR}/preassim_priorinf*
   if ( -e assimilate.csh )  ${REMOVE} ${RUN_DIR}/assimilate.csh

   echo "Listing contents of rundir before archiving at "`date`
   ls -l *.nc blown* dart_log* filter_* input.nml obs_seq* Output/inf_ic*
   mkdir -p ${OUTPUT_DIR}/${datea}/{Inflation_input,WRFIN,PRIORS,logs}

   set num_vars = $#increment_vars_a
   set extract_str = ''
   set i = 1
   while ( $i <= $num_vars )
      set extract_str = `echo ${extract_str}$increment_vars_a[$i],`
      @ i++
   end
   set extract_str = `echo ${extract_str}$increment_vars_a[$num_vars]`

   # Create an analysis increment file that has valid static data.
   # First, create the difference of a subset of variables
   # Second, create a netCDF file with just the static data
   # Third, append the static data onto the difference.
   ncdiff -F -O -v $extract_str postassim_mean.nc preassim_mean.nc analysis_increment.nc
   ncks -F -O -x -v ${extract_str} postassim_mean.nc static_data.nc
   ncks -A static_data.nc analysis_increment.nc

   # Move diagnostic and obs_seq.final data to storage directories

   foreach FILE ( postassim_mean.nc preassim_mean.nc postassim_sd.nc preassim_sd.nc \
                  obs_seq.final analysis_increment.nc output_mean.nc output_sd.nc )
      if ( -e $FILE && ! -z $FILE ) then
         ${MOVE} $FILE ${OUTPUT_DIR}/${datea}/.
         if ( ! $status == 0 ) then
            echo "failed moving ${RUN_DIR}/${FILE}"
            touch BOMBED
         endif
      else
         echo "${OUTPUT_DIR}/${FILE} does not exist and should."
         ls -l
         touch BOMBED
      endif
   end

   echo "past the analysis file moves"

   # Move inflation files to storage directories
   # The output inflation file is used as the input for the next cycle,
   # so rename the file 'on the fly'.
   cd ${RUN_DIR}   # TJH is this necessary?
   if ( $ADAPTIVE_INFLATION == 1 ) then
      set old_file = ( input_postinf_mean.nc  input_postinf_sd.nc  input_priorinf_mean.nc  input_priorinf_sd.nc )
      set new_file = ( output_postinf_mean.nc output_postinf_sd.nc output_priorinf_mean.nc output_priorinf_sd.nc )
      set i = 1
      set nfiles = $#new_file
      while ($i <= $nfiles)
         if ( -e ${new_file[$i]} && ! -z ${new_file[$i]} ) then
            ${MOVE} ${new_file[$i]} ${OUTPUT_DIR}/${datea}/Inflation_input/${old_file[$i]}
            if ( ! $status == 0 ) then
               echo "failed moving ${RUN_DIR}/Output/${FILE}"
               touch BOMBED
            endif
         endif
         @ i++
      end
      echo "past the inflation file moves"
   endif   # adaptive_inflation file moves

   # submit jobs to integrate ensemble members to next analysis time ...
   # BEFORE calculating the observation-space diagnostics for the existing cycle.

   echo "ready to integrate ensemble members"

   # Removing old start_member and done_member diagnostics
   if ( -e ${RUN_DIR}/start_member_1) then
      ${REMOVE} ${RUN_DIR}/start_member_*  \
                ${RUN_DIR}/done_member_*
   endif


   set n = 1
   while ( $n <= $NUM_ENS )

      if ( $SUPER_PLATFORM == 'LSF queuing system' ) then

         echo "2i\"                                                                  >! script.sed
         echo "#==================================================================\" >> script.sed
         echo "#BSUB -J assim_advance_${n}\"                                         >> script.sed
         echo "#BSUB -o assim_advance_${n}.%J.log\"                                  >> script.sed
         echo "#BSUB -P ${COMPUTER_CHARGE_ACCOUNT}\"                                 >> script.sed
         echo "#BSUB -W ${ADVANCE_TIME}\"                                            >> script.sed
         echo "#BSUB -q ${ADVANCE_QUEUE}\"                                           >> script.sed
         echo "#BSUB -n ${ADVANCE_CORES}\"                                           >> script.sed
         echo "#BSUB -x\"                                                            >> script.sed
         echo '#BSUB -R "span[ptile='"${ADVANCE_PTILE}"']"\'                         >> script.sed
         echo "#=================================================================="  >> script.sed
         echo 's%${1}%'"${datea}%g"                                                  >> script.sed
         echo 's%${2}%'"${n}%g"                                                      >> script.sed
         echo 's%${3}%'"${paramfile}%g"                                              >> script.sed

         sed -f script.sed ${SHELL_SCRIPTS_DIR}/assim_advance.csh >! assim_advance_mem${n}.csh
         if ( $?reservation ) then
            echo "MEMBER ${n} USING RESERVATION," `/contrib/lsf/get_my_rsvid`
            bsub -U `/contrib/lsf/get_my_rsvid` < assim_advance_mem${n}.csh
         else
            bsub < assim_advance_mem${n}.csh
         endif

      else if ( $SUPER_PLATFORM == 'derecho' ) then

         echo "2i\"                                                                             >! script.sed
         echo "#=================================================================\"             >> script.sed
         echo "#PBS -N assim_advance_${n}\"                                                     >> script.sed
         echo "#PBS -j oe\"                                                                     >> script.sed
         echo "#PBS -A ${COMPUTER_CHARGE_ACCOUNT}\"                                             >> script.sed
         echo "#PBS -l walltime=${ADVANCE_TIME}\"                                               >> script.sed
         echo "#PBS -q ${ADVANCE_QUEUE}\"                                                       >> script.sed
         echo "#PBS -l job_priority=${ADVANCE_PRIORITY}\"                                       >> script.sed
         echo "#PBS -m a\"                                                                      >> script.sed
         echo "#PBS -M ${EMAIL}\"                                                               >> script.sed
         echo "#PBS -k eod\"                                                                    >> script.sed
         echo "#PBS -l select=${ADVANCE_NODES}:ncpus=${ADVANCE_PROCS}:mpiprocs=${ADVANCE_MPI}\" >> script.sed
         echo "#================================================================="              >> script.sed
         echo 's%${1}%'"${datea}%g"                                                             >> script.sed
         echo 's%${2}%'"${n}%g"                                                                 >> script.sed
         echo 's%${3}%'"${paramfile}%g"                                                         >> script.sed

         sed -f script.sed ${SHELL_SCRIPTS_DIR}/assim_advance.csh >! assim_advance_mem${n}.csh
         qsub assim_advance_mem${n}.csh

      endif
      @ n++

   end

   #  Compute Diagnostic Quantities (in the background)

   if ( -e obs_diag.log ) ${REMOVE} obs_diag.log
   ${SHELL_SCRIPTS_DIR}/diagnostics_obs.csh $datea $paramfile >& ${RUN_DIR}/obs_diag.log &

   # ---------------------------------------------------------------------------
   #  check to see if all of the ensemble members have advanced
   # ---------------------------------------------------------------------------

   set advance_thresh = `echo $ADVANCE_TIME | cut -b3-4`
   @ advance_thresh = `expr $advance_thresh \+ 0` * 60 + `echo $ADVANCE_TIME | cut -b1-1` * 3600

   set n = 1
   while ( $n <= $NUM_ENS )

      set ensstring = `echo $n + 10000 | bc | cut -b2-5`
      set keep_trying = true

      while ( $keep_trying == 'true' )

         #  Wait for the script to start
         while ( ! -e ${RUN_DIR}/start_member_${n} )

            if ( $SUPER_PLATFORM == 'LSF queuing system' ) then

               if ( `bjobs -w | grep assim_advance_${n} | wc -l` == 0 ) then

                  echo "assim_advance_${n} is missing from the queue"
                  if ( $?reservation ) then
                     echo "MEMBER ${n} USING RESERVATION," `/contrib/lsf/get_my_rsvid`
                     bsub -U `/contrib/lsf/get_my_rsvid` < assim_advance_mem${n}.csh
                  else
                     bsub < assim_advance_mem${n}.csh
                  endif

               endif

            else if ( $SUPER_PLATFORM == 'derecho' ) then

               # Prevent double submission for member 1 only               
               if ( $n == 1) then
               sleep 5
               endif

               if ( `qstat -wa | grep assim_advance_${n} | wc -l` == 0 ) then

                  echo "assim_advance_${n} is missing from the queue"
                  qsub assim_advance_mem${n}.csh
               endif

            endif
            sleep 15

         end
         set start_time = `head -1 start_member_${n}`
         echo "Member $n has started.  Start time $start_time"

         #  Wait for the output file
         while ( 1 == 1 )

            set current_time = `date +%s`
            @ length_time = $current_time - $start_time

            if ( -e ${RUN_DIR}/done_member_${n} ) then

      	       #  If the output file already exists, move on
      	       set keep_trying = false
               break

            else if ( $length_time > $advance_thresh ) then

      	       #  Obviously, the job crashed.  Resubmit to queue
      	       ${REMOVE} start_member_${n}
               echo "didn't find the member done file"
               if ( $SUPER_PLATFORM == 'LSF queuing system' ) then

                  if ( $?reservation ) then
                     echo "MEMBER ${n} USING RESERVATION," `/contrib/lsf/get_my_rsvid`
      	             bsub -U `/contrib/lsf/get_my_rsvid` < assim_advance_mem${n}.csh
                  else
                     bsub < assim_advance_mem${n}.csh
      	          endif

               else if ( $SUPER_PLATFORM == 'derecho' ) then

                  qsub assim_advance_mem${n}.csh
                  sleep 5
               endif
      	       break

            endif
            sleep 15    # this might need to be longer, though I moved the done flag lower in the
                        # advance_model.csh to hopefully avoid the file moves below failing

         end

      end

      #  Move output data to correct location
      echo "moving ${n} ${ensstring}"
      ${MOVE} ${RUN_DIR}/assim_advance_${n}.o*              ${OUTPUT_DIR}/${datea}/logs/.
      ${MOVE} WRFOUT/wrf.out_${gdatef[1]}_${gdatef[2]}_${n} ${OUTPUT_DIR}/${datea}/logs/.
      ${MOVE} WRFIN/wrfinput_d01_${n}.gz                    ${OUTPUT_DIR}/${datea}/WRFIN/.
      ${MOVE} ${RUN_DIR}/prior_d01.${ensstring}             ${OUTPUT_DIR}/${datea}/PRIORS/.
      ${REMOVE} start_member_${n} done_member_${n} filter_restart_d01.${ensstring}
      if ( -e assim_advance_mem${n}.csh )  ${REMOVE} assim_advance_mem${n}.csh
      set pert = `cat ${RUN_DIR}/advance_temp${n}/mem${n}_pert_bank_num`
      echo "Member $n uses perturbation bank ensemble member $pert" >>  ${OUTPUT_DIR}/${datea}/pert_bank_members.txt

      @ n++

   end

   # ---------------------------------------------------------------------------
   # All ensemble members should have advanced by now.
   # ---------------------------------------------------------------------------

   if ( -e obs_prep.log ) ${REMOVE} obs_prep.log

      #  Clean everything up and finish

      #  Move DART-specific data to storage directory
      ${COPY} input.nml ${OUTPUT_DIR}/${datea}/.
      ${MOVE} ${RUN_DIR}/dart_log.out ${RUN_DIR}/dart_log.nml ${RUN_DIR}/*.log ${OUTPUT_DIR}/${datea}/logs/.

      #  Remove temporary files from both the run directory and old storage directories
      ${REMOVE} ${OUTPUT_DIR}/${datep}/wrfinput_d*_mean ${RUN_DIR}/wrfinput_d* ${RUN_DIR}/WRF

      #  Prep data for archive
      cd ${OUTPUT_DIR}/${datea}
      gzip -f wrfinput_d*_${gdate[1]}_${gdate[2]}_mean wrfinput_d*_${gdatef[1]}_${gdatef[2]}_mean wrfbdy_d*_mean
      tar -cvf retro.tar obs_seq.out wrfin*.gz wrfbdy_d*.gz
      tar -rvf dart_data.tar obs_seq.out obs_seq.final wrfinput_d*.gz wrfbdy_d*.gz \
                            Inflation_input/* logs/* *.dat input.nml
      ${REMOVE} wrfinput_d*_${gdate[1]}_${gdate[2]}_mean.gz wrfbdy_d*.gz
      gunzip -f wrfinput_d*_${gdatef[1]}_${gdatef[2]}_mean.gz

      cd $RUN_DIR
      ${MOVE} ${RUN_DIR}/assim*.o*            ${OUTPUT_DIR}/${datea}/logs/.
      ${MOVE} ${RUN_DIR}/*log                 ${OUTPUT_DIR}/${datea}/logs/.
      ${REMOVE} ${RUN_DIR}/input_priorinf_*
      ${REMOVE} ${RUN_DIR}/static_data*
      touch prev_cycle_done
      touch $RUN_DIR/cycle_finished_${datea}
      rm $RUN_DIR/cycle_started_${datea}

      # If doing a reanalysis, increment the time if not done.  Otherwise, let the script exit
      if ( $restore == 1 ) then
         if ( $datea == $datefnl) then
            echo "Reached the final date "
	    echo "Script exiting normally"
            exit 0
         endif
         set datea  = `echo $datea $ASSIM_INT_HOURS | ${DART_DIR}/models/wrf/work/advance_time`
      else
	 echo "Script exiting normally cycle ${datea}"
         exit 0
      endif
   end

exit 0

end

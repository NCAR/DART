#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# I am going to process one day per task and one month at a time via the job
# array syntax.
#
#BXXX -J qscat[305-334]   # nov 2007
#BXXX -J qscat[335-365]   # dec 2007
#BXXX -J qscat[1-31]      # jan 2008
#BXXX -J qscat[32-60]     # feb 2008
#BXXX -J qscat[61-92]     # mar 2008
#
#BSUB -J qscat[335-365]
#BSUB -n 1
#BSUB -q standby
#BSUB -W 0:10
#BSUB -o qscat.%J.%I
#BSUB -N -u ${USER}@ucar.edu
#BSUB -m  "cr0128en cr0129en cr0130en cr0131en cr0132en cr0133en cr0134en cr0135en cr0136en cr0137en cr0138en cr0139en cr0140en cr0141en cr0202en cr0201en"
#
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------
#----------------------------------------------------------------------

set month = 12

if ($?LSB_HOSTS) then

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_OUTPUTFILE:ar
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST
   setenv TASKID      $LSB_JOBINDEX
   setenv NTASKS      $LSB_JOBINDEX_END
   setenv STEP        $LSB_JOBINDEX_STEP

   set DARTHOME = /fs/image/home/thoar/DART/observations
   set DATADIR = /ptmp/thoar/QuikSCAT_L2B

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     QuikSCAT
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $host
   setenv TASKID      305
   setenv NTASKS      1
   setenv STEP        1

   set DARTHOME = /fs/image/home/thoar/DART/observations
   set DATADIR = /ptmp/thoar/QuikSCAT_L2B

endif

echo ""
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) running on  host $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) job $TASKID of $NTASKS (by $STEP) started at "`date`
echo ""

#----------------------------------------------------------------------
#----------------------------------------------------------------------

@ y = 2007
@ d = 305
@ dmax = 313

while ( $d <= $dmax )

   cd /ptmp/thoar/QuikSCAT_L2B/${y}/${d}

   \rm -f *obs_seq_out

   gunzip -q QS_S2B*.gz

   foreach FILE ( QS_S2B* )

      echo "&obs_kind_nml"                             >! input.nml
      echo "   /"                                      >> input.nml
      echo "&location_nml"                             >> input.nml
      echo "   /"                                      >> input.nml
      echo "&utilities_nml"                            >> input.nml
      echo "   /"                                      >> input.nml
      echo "&obs_sequence_nml"                         >> input.nml
      echo "   write_binary_obs_sequence = .false."    >> input.nml
      echo "   /"                                      >> input.nml
      echo " "                                         >> input.nml
      echo "&convert_L2b_nml"                          >> input.nml
      echo "datadir   = '.',"                          >> input.nml
      echo "outputdir = '.',"                          >> input.nml
      echo "l2b_file = '"$FILE"',"                     >> input.nml
      echo "lon1 =   0.0, lon2 = 360.0,"               >> input.nml
      echo "lat1 = -90.0, lat2 =  90.0"                >> input.nml
      echo "   /"                                      >> input.nml
      echo " "                                         >> input.nml

      echo "      $FILE"

      ${DARTHOME}/quikscat/work/convert_L2b > /dev/null || exit 1
   
   end

   #----------------------------------------------------------------------
   # Concatenate all the observation sequence files for each orbit 
   # into one observation sequence file for a day. 
   #
   # Create namelist for obs_sequence_tool
   #----------------------------------------------------------------------

   set STRING = "1,$ s# #', '#g"
   set filenames = `ls *obs_seq_out`
   set numorbits = $#filenames

   if ( $numorbits > 0 ) then

      echo $filenames >! filenames_file

      set filenamestring = `sed -e "$STRING" filenames_file` 

      echo " &obs_sequence_tool_nml"                >> input.nml
      echo " num_input_files = ${numorbits},"       >> input.nml
      echo " filename_seq = '"${filenamestring}"'," >> input.nml
      echo " filename_out = 'obs_seq.processed',"   >> input.nml
      echo " first_obs_days = -1,"                  >> input.nml
      echo " first_obs_seconds = -1,"               >> input.nml
      echo " last_obs_days = -1,"                   >> input.nml
      echo " last_obs_seconds = -1,"                >> input.nml
      echo " obs_types = '',"                       >> input.nml
      echo " keep_types = .false.,"                 >> input.nml
      echo " print_only = .false.,"                 >> input.nml
      echo " min_lat = -90.0,"                      >> input.nml
      echo " max_lat = 90.0,"                       >> input.nml
      echo " min_lon = 0.0,"                        >> input.nml
      echo " max_lon = 360.0"                       >> input.nml
      echo "    /"                                  >> input.nml

      ${DARTHOME}/utilities/threed_sphere/obs_sequence_tool || exit 2

      mv obs_seq.processed ${DATADIR}/${y}11/qscatL2B_${y}_${d}_obs_seq.out

   endif

   @ d = $d + 1

   \rm -rf input.nml filenames_file dart_log.out dart_log.nml

end

echo "${JOBNAME} ($JOBID) job $TASKID finished at "`date`

exit 0



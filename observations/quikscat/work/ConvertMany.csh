#!/bin/csh
#
#BSUB -J qscat
#BSUB -n 1
#BSUB -q standby
#BSUB -W 0:05
#BSUB -o out.%J.%I
#BSUB -e err.%J.%I
#BSUB -N -u ${USER}@ucar.edu

#----------------------------------------------------------------------
#----------------------------------------------------------------------

set DARTHOME = /fs/image/home/thoar/DART/observations
set DATADIR = /gpfs/ptmp/dart/Obs_sets/QuikSCAT_24_ascii

#----------------------------------------------------------------------
#----------------------------------------------------------------------

@ y = 2007
@ n = 305
@ nmax = 313

while ( $n <= $nmax )

   echo "$y $n"

   cd /ptmp/thoar/QuikSCAT_L2B/${y}/${n}

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

      mv obs_seq.processed ${DATADIR}/${y}11/qscatL2B_${y}_${n}_obs_seq.out

   endif

   @ n = $n + 1

   \rm -rf input.nml filenames_file dart_log.out dart_log.nml

end


#!/bin/ksh 

# test suite for merge_obs_sequence utility

# set up constant part at beginning
cat > input.nml.head <<EOF
&obs_sequence_nml
   write_binary_obs_sequence = .false.  /

&obs_kind_nml
   assimilate_these_obs_types = 'RAW_STATE_VARIABLE'  /

&utilities_nml
   TERMLEVEL = 1,
   module_details = .false.,
   logfilename = 'dart_log.out'  /


EOF

# we will concatinate onto the end of this file
rm -fr test_merge.log

# should have a second arg - expects it to fail or not

# set up functions which will be used for the rest of the script
# shortone cats the entire result into the log, longone heads the
# first 30 and last 30 lines
function shortone
{
   rm -f input.nml
   cat input.nml.head input.fragment > input.nml
   echo "" >> test_merge.log
   echo "-------------------------" >> test_merge.log
   echo "-------------------------" >> test_merge.log
   echo " test $1 " >> test_merge.log
   echo "-------------------------" >> test_merge.log
   echo "" >> test_merge.log
   cat input.fragment >> test_merge.log
   echo "" >> test_merge.log
   echo "-------------------------" >> test_merge.log
   rm -f obs_seq.merged
   ./merge_obs_seq 2>&1 >> test_merge.log
   echo "-------------------------" >> test_merge.log
   if [[ -f obs_seq.merged ]]; then
     cat obs_seq.merged >> test_merge.log
   else
     echo obs_seq.merged not found  >> test_merge.log
   fi
}

function longone
{
   rm -f input.nml
   cat input.nml.head input.fragment > input.nml
   echo "" >> test_merge.log
   echo "-------------------------" >> test_merge.log
   echo "-------------------------" >> test_merge.log
   echo " test $1 " >> test_merge.log
   echo "-------------------------" >> test_merge.log
   echo "" >> test_merge.log
   cat input.fragment >> test_merge.log
   echo "" >> test_merge.log
   echo "-------------------------" >> test_merge.log
   rm -f obs_seq.merged
   ./merge_obs_seq 2>&1 >> test_merge.log
   echo "-------------------------" >> test_merge.log
   if [[ -f obs_seq.merged ]]; then
     head -30 obs_seq.merged >> test_merge.log
     echo '-------------------------' >> test_merge.log
     tail -30 obs_seq.merged >> test_merge.log
     echo '-------------------------' >> test_merge.log
   else
     echo obs_seq.merged not found  >> test_merge.log
   fi
}

###############################
cat - > input.fragment <<EOF

&merge_obs_seq_nml
   num_input_files = 2,
   filename_seq    = 'obs_seq.A', 'obs_seq.A',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = -1,
   first_obs_seconds        = -1,
   last_obs_days            = -1,
   last_obs_seconds         = -1,
  /

EOF
shortone with_self works


###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 2,
   filename_seq    = 'obs_seq.A', 'obs_seq.A',
   filename_out    = 'obs_seq.merged',
  /
EOF

shortone with_self_no_times works


###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 1,
   filename_seq    = 'obs_seq.A',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = -1,
   first_obs_seconds        = -1,
   last_obs_days            = -1,
   last_obs_seconds         = -1,
  /
EOF

shortone self_only works


###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 2,
   filename_seq    = 'obs_seq.A', 'obs_seq.B',
   filename_out    = 'obs_seq.merged',
  /
EOF

shortone should_fail fails

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 2,
   filename_seq    = 'obs_seq.A', 'obs_seq.C',
   filename_out    = 'obs_seq.merged',
  /
EOF

shortone should_fail fails

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 2,
   filename_seq    = 'obs_seq.A', 'obs_seq.D',
   filename_out    = 'obs_seq.merged',
  /
EOF

shortone should_fail fails

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 2,
   filename_seq    = 'obs_seq.A', 'obs_seq.E',
   filename_out    = 'obs_seq.merged',
  /
EOF

shortone should_fail fails

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 2,
   filename_seq    = 'obs_seq.A', 'obs_seq.F',
   filename_out    = 'obs_seq.merged',
  /
EOF

shortone should_fail fails

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 2,
   filename_seq    = 'obs_seq.A', 'obs_seq.G',
   filename_out    = 'obs_seq.merged',
  /
EOF

shortone should_work works


###############################
###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 3,
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
  /
EOF

longone case_0 works

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 3,
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 0,
   last_obs_days            = 0,
   last_obs_seconds         = 4000,
  /
EOF

longone case_1 works

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 3,
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 7000,
   last_obs_days            = 0,
   last_obs_seconds         = 11000,
  /
EOF

longone case_2 works

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 3,
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 20000,
   last_obs_days            = 0,
   last_obs_seconds         = 25000,
  /
EOF

longone case_3 works


###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 3,
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 15000,
   last_obs_days            = 0,
   last_obs_seconds         = 16000,
  /
EOF

longone case_4 fails

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 3,
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 0,
   last_obs_days            = 0,
   last_obs_seconds         = 2000,
  /
EOF

longone case_5 fails

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 3,
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 25000,
   last_obs_days            = 0,
   last_obs_seconds         = 28000,
  /
EOF

longone case_6 fails


###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 4,
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 0,
   last_obs_days            = 0,
   last_obs_seconds         = 10000,
  /
EOF

longone case_7 works


###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 4,
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 13000,
   last_obs_days            = 0,
   last_obs_seconds         = 25000,
  /
EOF

longone case_8 works


###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 4,
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 4000,
   last_obs_days            = 0,
   last_obs_seconds         = 7000,
  /
EOF

longone case_9 works

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 4,
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 7200,
   last_obs_days            = 0,
   last_obs_seconds         = 9000,
  /
EOF

longone case_9a works

###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 4,
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 10000,
   last_obs_days            = 0,
   last_obs_seconds         = 11000,
  /
EOF

longone case_9b works


###############################
cat > input.fragment <<EOF
&merge_obs_seq_nml
   num_input_files = 4,
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 10000,
   last_obs_days            = 0,
   last_obs_seconds         = 18000,
  /
EOF

longone case_10 works



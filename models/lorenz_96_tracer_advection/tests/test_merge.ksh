#!/bin/ksh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# test suite for merge_obs_sequence utility

# set up constant part at beginning
cat > input.nml.head <<EOF
&obs_sequence_nml
   write_binary_obs_sequence = .false.  
   /

&obs_kind_nml
   /

&utilities_nml
   TERMLEVEL = 1,
   module_details = .false.,
   /

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
   echo " test $1 ; expect $2" >> test_merge.log
   echo "" >> test_merge.log
   cat input.fragment >> test_merge.log
   echo "" >> test_merge.log
   rm -f obs_seq.merged
   ./obs_sequence_tool 2>&1 >> test_merge.log
   if [[ -f obs_seq.merged ]]; then
     if [[ $2 == 'fails' ]]; then
       echo TESTERROR: test succeeded when it should have failed >> test_merge.log
     fi
     #cat obs_seq.merged >> test_merge.log
   else
     if [[ $2 == 'works' ]]; then
       echo TESTERROR: test failed when it should have worked >> test_merge.log
     fi
     #echo obs_seq.merged not found  >> test_merge.log
   fi
}

function longone
{
   rm -f input.nml
   cat input.nml.head input.fragment > input.nml
   echo "" >> test_merge.log
   echo "-------------------------" >> test_merge.log
   echo " test $1 ; expect $2" >> test_merge.log
   echo "" >> test_merge.log
   cat input.fragment >> test_merge.log
   echo "" >> test_merge.log
   rm -f obs_seq.merged
   ./obs_sequence_tool 2>&1 >> test_merge.log
   if [[ -f obs_seq.merged ]]; then
     if [[ $2 == 'fails' ]]; then
       echo TESTERROR: test succeeded when it should have failed >> test_merge.log
     fi
     #head -30 obs_seq.merged >> test_merge.log
     #echo '-------------------------' >> test_merge.log
     #tail -30 obs_seq.merged >> test_merge.log
     #echo '-------------------------' >> test_merge.log
   else
     if [[ $2 == 'works' ]]; then
       echo TESTERROR: test failed when it should have worked >> test_merge.log
     fi
     #echo obs_seq.merged not found  >> test_merge.log
   fi
}

###############################
cat - > input.fragment <<EOF

&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.A', 'obs_seq.A',
   filename_out    = 'obs_seq.merged',
   gregorian_cal   = .false.
  /

EOF
shortone with_self works


###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.A', 'obs_seq.A',
   filename_out    = 'obs_seq.merged',
   gregorian_cal   = .false.
  /
EOF

shortone with_self_no_times works


###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.A',
   filename_out    = 'obs_seq.merged',
   gregorian_cal     = .false.
  /
EOF

shortone self_only works


###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.A', 'obs_seq.B',
   filename_out    = 'obs_seq.merged',
   gregorian_cal   = .false.
  /
EOF

shortone should_fail1 fails

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.A', 'obs_seq.C',
   filename_out    = 'obs_seq.merged',
   gregorian_cal   = .false.
  /
EOF

shortone should_fail2 fails

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.A', 'obs_seq.D',
   filename_out    = 'obs_seq.merged',
   gregorian_cal   = .false.
  /
EOF

shortone should_fail3 fails

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.A', 'obs_seq.E',
   filename_out    = 'obs_seq.merged',
   gregorian_cal   = .false.
  /
EOF

shortone should_work1 works

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.A', 'obs_seq.F',
   filename_out    = 'obs_seq.merged',
   gregorian_cal   = .false.
  /
EOF

shortone should_fail5 fails

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.A', 'obs_seq.G',
   filename_out    = 'obs_seq.merged',
   gregorian_cal   = .false.
  /
EOF

shortone should_work2 works


###############################
###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   gregorian_cal   = .false.
  /
EOF

longone case_0 works

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 0,
   last_obs_days            = 0,
   last_obs_seconds         = 4000,
   gregorian_cal   = .false.
  /
EOF

longone case_1 works

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 7000,
   last_obs_days            = 0,
   last_obs_seconds         = 11000,
   gregorian_cal   = .false.
  /
EOF

longone case_2 works

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 20000,
   last_obs_days            = 0,
   last_obs_seconds         = 25000,
   gregorian_cal   = .false.
  /
EOF

longone case_3 works


###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 15000,
   last_obs_days            = 0,
   last_obs_seconds         = 16000,
   gregorian_cal   = .false.
  /
EOF

longone case_4 fails

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 0,
   last_obs_days            = 0,
   last_obs_seconds         = 2000,
   gregorian_cal   = .false.
  /
EOF

longone case_5 fails

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.L', 'obs_seq.M', 'obs_seq.N',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 25000,
   last_obs_days            = 0,
   last_obs_seconds         = 28000,
   gregorian_cal   = .false.
  /
EOF

longone case_6 fails


###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 0,
   last_obs_days            = 0,
   last_obs_seconds         = 10000,
   gregorian_cal   = .false.
  /
EOF

longone case_7 works


###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 13000,
   last_obs_days            = 0,
   last_obs_seconds         = 25000,
   gregorian_cal   = .false.
  /
EOF

longone case_8 works


###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 4000,
   last_obs_days            = 0,
   last_obs_seconds         = 7000,
   gregorian_cal   = .false.
  /
EOF

longone case_9 works

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 7200,
   last_obs_days            = 0,
   last_obs_seconds         = 9000,
   gregorian_cal   = .false.
  /
EOF

longone case_9a works

###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 10000,
   last_obs_days            = 0,
   last_obs_seconds         = 11000,
   gregorian_cal   = .false.
  /
EOF

longone case_9b works


###############################
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq    = 'obs_seq.O', 'obs_seq.P', 'obs_seq.Q', 'obs_seq.R',
   filename_out    = 'obs_seq.merged',
   first_obs_days           = 0,
   first_obs_seconds        = 10000,
   last_obs_days            = 0,
   last_obs_seconds         = 18000,
   gregorian_cal   = .false.
  /
EOF

longone case_10 works

###############################
rm -f flist
echo obs_seq.O obs_seq.P obs_seq.Q obs_seq.R | xargs -n 1 echo >> flist
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq      = ''
   filename_seq_list = 'flist'
   filename_out      = 'obs_seq.merged',
   gregorian_cal     = .false.
  /
EOF

longone case_11 works

###############################
rm -f flist
echo obs_seq.O obs_seq.P obs_seq.Q obs_seq.R | xargs -n 1 echo >> flist
cat > input.fragment <<EOF
&obs_sequence_tool_nml
   filename_seq      = 'obs_seq.O'
   filename_seq_list = 'flist'
   filename_out      = 'obs_seq.merged',
   gregorian_cal     = .false.
  /
EOF

longone case_12 fails

###############################
###############################

n=`fgrep TESTERROR: test_merge.log | wc -l`
echo found $n bad tests
if [[ $n -gt 0 ]]; then
  echo check output in log file:  test_merge.log
  echo bad tests are marked with \'TESTERROR:\'
fi

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


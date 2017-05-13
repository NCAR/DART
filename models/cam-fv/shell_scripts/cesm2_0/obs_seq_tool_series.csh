#!/bin/tcsh

# Script to process a series of obs_seq.final files,
# to change any of the properties available to obs_sequence_tool.
#
# 1) copy input.nml to the directory containing the subdirectory of obs_seq files.
# 2) edit input.nml:obs_sequence_tool_nml to modify as needed.
#       filename_seq       = 'obs_seq.in',
#       filename_out       = 'obs_seq.out',
#
# 3) edit this script as needed.
#

set tool = ~/DART/Trunk/models/cam/work/obs_sequence_tool

set in_dir  = Obs_seqs
set out_dir = Obs_seqs_nomembers

if (! -d $in_dir) then
   echo "Missing $in_dir"
   exit 5
endif

if (-d $out_dir) then
   echo "$out_dir already exists; move it or pick another name"
   exit 10
else
   mkdir $out_dir
endif

foreach f (`ls $in_dir`)
   ln -s $in_dir/$f obs_seq.in
   $tool
   if (-f obs_seq.out) then
      mv obs_seq.out $out_dir/$f
      rm obs_seq.in
   else
      echo "no obs_seq.out from $f.  Aborting"
      exit 20
   endif
end


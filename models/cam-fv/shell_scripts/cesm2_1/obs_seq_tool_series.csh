#!/bin/tcsh

# Script to process a series of obs_seq.final files,
# to change any of the properties available to obs_sequence_tool.
# > > > In this case, remove AIRS obs.
#
# 1) copy input.nml to the directory containing the subdirectory of obs_seq files.
# 2) edit input.nml:obs_sequence_tool_nml to modify as needed.
#       filename_seq       = 'obs_seq.in',
#       filename_out       = 'obs_seq.out',
#
# 3) edit this script as needed.
#


# should this script reference data_scripts.csh for these values?
if ($#argv == 0) then
   set case_base = f.e21.FHIST_BGC.f09_025.CAM6assim
   set case = 004
else
   set case_base = $1
   set case = $2
endif

set tool = ~/DART/reanalysis/models/cam-fv/work_casper/obs_sequence_tool

# set in_dir  = Obs_seqs
# set out_dir = Obs_seqs_nomembers
set in_dir  = /glade/scratch/raeder/${case_base}.${case}/archive/esp/hist
set out_dir = Obs_${case}_noAIRS

if (! -d $in_dir) then
   echo "obs_seq_tool_series.csh: Missing $in_dir"
   exit 5
endif

if (-d $out_dir) then
   echo "obs_seq_tool_series.csh: $out_dir already exists; move it or pick another name"
   exit 10
else
   mkdir $out_dir
endif

# f.e21.FHIST_BGC.f09_025.CAM6assim.002.dart.e.cam_obs_seq_final.2017-01-01-21600
# foreach f (`ls $in_dir/*obs_seq_final*`)
foreach f (`cat obs_${case}.list`)
   ln -sf $in_dir/$f obs_seq.in
   echo "obs_seq_tool_series.csh tool ="
   ls -l $tool

   $tool >& obs_seq_tool.eo
   if (${status} != 0) then
      echo "Failed with status $status"
      exit 20
   endif

   if (-f obs_seq.out) then
      mv obs_seq.out $out_dir/$f
      rm obs_seq.in
   else
      echo "obs_seq_tool_series.csh: no obs_seq.out from $f.  Aborting"
      exit 20
   endif
end

# The /* causes the directory name to be included in the output.
# Then the files can be found from $diagnostics_dir in diags_batch.csh.
ls -1 $out_dir/* >! ${out_dir}.list

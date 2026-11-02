#!/bin/tcsh

# After the core is working, expand to do all components.
# Move the existing $comp/work/input.nml to input.nml.not_cpld
# set cesm = $1
# set parts = (`echo $cesm | sed -e 's#\.# #g'`)
# set comps = (`echo $parts[3] | sed -e 's#+# #g'`)
# 
# foreach comp ($comps)
#    ... the stuff below, with 
#        input.nml.cam-fv+clm -> $cesm, 
#           The extractions from the cesm input.nml can be preserved for use by all components
#        cam-fv -> $comp
# end
# 
# Doing them separately puts them in a predictable order.
# for testing all files are in CESM/work/.  Normally the tor file will be in its work dir.
set echo
# Extract from the CESM input.nml
foreach nml (preprocess obs_kind)
   grep -n -m 1 -A 100 -e "&${nml}_nml"  input.nml.cam-fv+clm \
     | sed -e "s#[:-]# #" >! batch
   set line = (`head -n 1 batch`)
   @ chunk_start = $line[1]
   set line = (`grep -m 1 -e ' /' batch `)
   @ chunk_end = $line[1] 

   # Keep only the namelist and replace the line numbers with spaces
   sed -n $chunk_start,${chunk_end}p input.nml.cam-fv+clm \
       | sed  -e "s#^[0-9][0-9 ][ ]?#  #" >! ${nml}.chunk
end
# The chunks have exactly what I want from cam-fv+clm

# preprocess tor
grep -n -A 100 -e "&preprocess_nml"  input.nml.cam-fv | sed -e "s#[:-]# #" >! batch.prep
set line = (`head -n 1 batch.prep`)
@ prep_start = $line[1]

# obs_kind tor
grep -n -A 100 -e "&obs_kind_nml"  input.nml.cam-fv | sed -e "s#[:-]# #" >! batch.obs
set line = (`head -n 1 batch.obs`)
@ obs_start = $line[1]

# Might need to rearrange preproc and obs_kind order
if ($prep_start < $obs_start) then
   @ fill_end = $prep_start - 1
   echo fill1
   head -n $fill_end input.nml.cam-fv >! fill1

   echo fill2
   set line = (`grep -m 1 ' /' batch.prep`)
   @ prep_end = $line[1]
   @ fill_start = $prep_end + 1
   @ fill_end = $obs_start - 1
   # -n prevents duplicate lines in the range selected
   #    and the printing (once) of all non-selected lines.
   sed -n $fill_start,${fill_end}p input.nml.cam-fv >! fill2

   echo fill3
   set line = (`grep -m 1 ' /' batch.obs`)
   @ obs_end = $line[1]
   @ fill_start = $obs_end + 1
   set line = (`wc -l input.nml.cam-fv`)
   @ fill_end = $line[1]
   sed -n $fill_start,${fill_end}p input.nml.cam-fv >! fill3

else
   @ fill_end = $obs_start - 1
   head -n $fill_end input.nml.cam-fv >! fill1

   set line = (`grep -m 1 ' /' batch.obs`)
   @ obs_end = $line[1]
   @ fill_start = $obs_end + 1
   @ fill_end = $prep_start - 1
   sed -n $fill_start,${fill_end}p input.nml.cam-fv >! fill2

   set line = (`grep -m 1 ' /' batch.prep`)
   @ prep_end = $line[1]
   @ fill_start = $prep_start + $prep_end 
   set line = `wc -l input.nml.cam-fv`	
   @ fill_end = $line[1]
   sed -n $fill_start,${fill_end}p input.nml.cam-fv >! fill3

endif

cat fill1 preprocess.chunk fill2 obs_kind.chunk fill3 >! input.nml

# rm fill* batch* 


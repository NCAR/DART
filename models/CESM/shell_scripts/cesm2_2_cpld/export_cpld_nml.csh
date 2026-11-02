#!/bin/tcsh

# The core is working (export_cpld_nml_cam-fv.csh), so expand to do all components.

# This script inherits useful parameters (env vars) from setup_assimilation
# and DART_params.

# for testing all files are in CESM/work/.  Normally the tor file will be in its work dir.
set echo
# Extract from the CESM input.nml
# Doing them separately puts them in a predictable order.
foreach nml (preprocess obs_kind)
   grep -n -m 1 -A 100 -e "&${nml}_nml"  ${shared_input_nml} \
     | sed -e "s#[:-]# #" >! batch
   set line = (`head -n 1 batch`)
   @ chunk_start = $line[1]
   set line = (`grep -m 1 -e ' /' batch `)
   @ chunk_end = $line[1] 

   # Keep only the namelist and replace the line numbers with spaces
   sed -n $chunk_start,${chunk_end}p ${shared_input_nml} \
       | sed  -e "s#^[0-9][0-9 ][ ]?#  #" >! ${nml}.chunk
end
# The chunks have exactly what I want from cam-fv+clm

# Break the input.nml of each assimilator into sections,
# so that preprocess_nml and obs_kind_nml can be replaced:
#   proprocess_nml
#   obs_kind_nml
#   all the text before, between, and after.
# This preserves the structure of the assimilator's input.nml
# so that it can be compared easily to other versions.
# This is run from CESM/work (directories were changed in setup_assimilation)
foreach tor ($assimilators)
   set tor_input = ../../${tor}/work/input.nml
   if (! -f $tor_input) then
      echo "ERROR: $0 is missing $tor_input"
      exit 77
   endif
   $COPY $tor_input ${tor_input}.not_cpld
   $MOVE $tor_input input.nml.cp

   # The preprocess of the tor
   grep -n -A 100 -e "&preprocess_nml"  input.nml.cp \
        | sed -e "s#[:-]# #" >! batch.prep
   set line = (`head -n 1 batch.prep`)
   @ prep_start = $line[1]
   
   # The obs_kind of the tor
   grep -n -A 100 -e "&obs_kind_nml"  input.nml.cp \
        | sed -e "s#[:-]# #" >! batch.obs
   set line = (`head -n 1 batch.obs`)
   @ obs_start = $line[1]
   
   # Might need to rearrange preproc and obs_kind order
   if ($prep_start < $obs_start) then
      @ fill_end = $prep_start - 1
      # echo fill1
      head -n $fill_end input.nml.cp >! fill1
   
      # echo fill2
      set line = (`grep -m 1 ' /' batch.prep`)
      @ prep_end = $line[1]
      @ fill_start = $prep_end + 1
      @ fill_end = $obs_start - 1
      # -n prevents duplicate lines in the range selected
      #    and the printing (once) of all non-selected lines.
      sed -n $fill_start,${fill_end}p input.nml.cp >! fill2
   
      # echo fill3
      set line = (`grep -m 1 ' /' batch.obs`)
      @ obs_end = $line[1]
      @ fill_start = $obs_end + 1
      set line = (`wc -l input.nml.cp`)
      @ fill_end = $line[1]
      sed -n $fill_start,${fill_end}p input.nml.cp >! fill3
   
   else
      @ fill_end = $obs_start - 1
      head -n $fill_end input.nml.cp >! fill1
   
      set line = (`grep -m 1 ' /' batch.obs`)
      @ obs_end = $line[1]
      @ fill_start = $obs_end + 1
      @ fill_end = $prep_start - 1
      sed -n $fill_start,${fill_end}p input.nml.cp >! fill2
   
      set line = (`grep -m 1 ' /' batch.prep`)
      @ prep_end = $line[1]
      @ fill_start = $prep_start + $prep_end 
      set line = (`wc -l input.nml.cp`)
      @ fill_end = $line[1]
      sed -n $fill_start,${fill_end}p input.nml.cp >! fill3
   
   endif
end

cat fill1 preprocess.chunk fill2 obs_kind.chunk fill3 >! $tor_input

# rm fill* batch* 


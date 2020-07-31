#!/bin/csh 

# generate list dynamically, sed out the string ...

set MYSTRING = :observations/forward_operators/obs_def_mod.f90

foreach FILE ( `grep obs_def_mod.f90 path* | sed -e "s#${MYSTRING}##"` )
  echo "Adding rttov module to $FILE"  
  echo "observations/forward_operators/rttov_interface_mod.f90" >> $FILE
end

echo ""
echo "Make sure that you have added 'observations/forward_operators/obs_def_rttov_mod.f90'"
echo "to the input.nml  preprocess_nml:input_files variable to support the use of the"
echo "RTTOV forward operators. Naturally, you must run 'preprocess' before you run any"
echo "of the individual mkmf_xxxx files. quickbuild.csh automatically does this."
echo ""


#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

if ( $#argv > 0 ) then
   set n     = ${1}   # pass in the ensemble member number
   set datep = ${2}   # needed for correct path to file
   set dn    = ${3}
   set paramfile = ${4}
else # values come from environment variables   #TJH If these are not set ....
   set n     = $mem_num
   set datep = $date
   set dn    = $domain
   set paramfile = $paramf
endif
source $paramfile

echo "prep_ic.csh using n=$n datep=$datep dn=$dn paramfile=$paramf"

if ( $dn == 1 ) then

   set num_vars = $#cycle_vars_a    # defined in paramfile
   set cycle_str = ''   # these are variables we want to cycle
   set i = 1
   while ( $i < $num_vars )
      set cycle_str = `echo ${cycle_str}$cycle_vars_a[$i],`
      @ i ++
   end
   set cycle_str = `echo ${cycle_str}$cycle_vars_a[$num_vars]`
   echo ${cycle_str}

else   # larger domain numbers use a different list of cycled variables (includes radar)

   set num_vars = $#cycle_vars_b    # defined in paramfile
   set cycle_str = ''   # these are variables we want to cycle
   set i = 1
   while ( $i < $num_vars )
      set cycle_str = `echo ${cycle_str}$cycle_vars_b[$i],`
      @ i ++
   end
   set cycle_str = `echo ${cycle_str}$cycle_vars_b[$num_vars]`
   echo ${cycle_str}

endif

set ensstring = `printf %04d $n`
set dchar     = `printf %02d $dn`

ncks -A -v ${cycle_str} \
          ${OUTPUT_DIR}/${datep}/PRIORS/prior_d${dchar}.${ensstring} \
          ${RUN_DIR}/advance_temp${n}/wrfinput_d${dchar}

touch ${RUN_DIR}/ic_d${dchar}_${n}_ready

exit 0


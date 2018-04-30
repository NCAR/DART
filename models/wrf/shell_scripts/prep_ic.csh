#!/bin/csh
########################################################################
  set paramfile = /glade2/scratch2/USERNAME/WORK_DIR/scripts/param.csh   # set this appropriately #%%%#
  source $paramfile

  if ( $#argv > 0 ) then 
   set n = ${1}      # pass in the ensemble member number
   set datep = ${2}  # needed for correct path to file
   set dn    = ${3}
  else
   set n = $mem_num      # pass in the ensemble member number
   set datep = $date  # needed for correct path to file
   set dn    = $domain
  endif

  echo "$n  $datep  $dn"

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

  set ensstring = `echo $n + 10000 | bc | cut -b2-5`
  set dchar     = `echo $dn + 100 | bc | cut -b2-3`

  ncks -A -v ${cycle_str} ${OUTPUT_DIR}/${datep}/PRIORS/prior_d${dchar}.${ensstring} ${RUN_DIR}/advance_temp${n}/wrfinput_d${dchar}

  touch ${RUN_DIR}/ic_d${dchar}_${n}_ready

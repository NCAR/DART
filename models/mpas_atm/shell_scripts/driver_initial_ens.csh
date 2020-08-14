#!/bin/csh
####################################################################################
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#   driver_initial_ens.csh - driver shell script that generates an
#                            initial MPAS ensemble from an ensemble of
#                            grib files.  Job will submit individual
#                            members to a queuing system (or run them
#                            sequentially).
####################################################################################

if ( $#argv >= 1 ) then
   set fn_param = ${1}
else
   set fn_param = `pwd`/setup_params.csh
endif

if(! -e $fn_param ) then
   echo $fn_param does not exist. Cannot proceed.
   exit
endif

source $fn_param

cd $RUN_DIR

foreach fn ( $INIT_GRIB_FILE_LIST )
  if ( ! -r ${fn} ) then
    echo ABORT\: could not find required readable dependency ${fn}
    exit 1
  endif 
end

foreach fn ( advance_time )
   if ( ! -x $fn ) then
      echo ${COPY} ${DART_DIR}/${fn} .
           ${COPY} ${DART_DIR}/${fn} .
      if ( ! $status == 0 ) then
         echo ABORT\: We cannot find required executable dependency $fn.
         exit
      endif
   endif
end

foreach fn ( init_mpas_grib.csh )
  if ( ! -r $fn ) then
     echo ${COPY} ${DART_DIR}/../shell_scripts/${fn} .
          ${COPY} ${DART_DIR}/../shell_scripts/${fn} .
     if ( ! $status == 0 ) then
        echo ABORT\: We cannot find required script $fn.
        exit
     endif
  endif
end

if ( ! -d MPAS_RUN ) then

   if ( ! -d $MPAS_DIR ) then
      echo $MPAS_DIR does not exist. Stop.
      exit
   endif
   ${LINK} $MPAS_DIR MPAS_RUN

endif

#  Check to see if MPAS and DART namelists exist.  If not, copy them from model source
foreach fn ( ${NML_MPAS} ${NML_INIT} )
  if ( ! -r ${fn} ) then
    ${COPY} MPAS_RUN/${fn} .
  endif
end

foreach fn ( ${STREAM_ATM} ${STREAM_INIT} )
  if ( ! -r ${fn} || -z $fn ) then
    ${COPY} MPAS_RUN/${fn} .
  endif
end

if ( ! -r ${NML_DART} ) then
  ${COPY} ${DART_DIR}/${NML_DART} .
endif

if ( ! -x ungrib.exe ) then
  ${COPY} ${WPS_DIR}/ungrib.exe .
endif

if ( ! -r ${NML_WPS} ) then
  ${COPY} ${WPS_DIR}/${NML_WPS} .
endif

if ( ! -r Vtable ) then
  ${COPY} ${WPS_DIR}/${VTABLE} Vtable
endif

@ ndecomp = $MODEL_NODES * $N_PROCS
set fgraph = ${MPAS_GRID}.graph.info.part.${ndecomp}
if ( ! -e ${fgraph} ) then
   ${LINK} ${GRID_DIR}/${fgraph} ${fgraph}
   if(! -e ${fgraph}) then
      echo "Cannot find ${fgraph} for n_mpas * n_proc (= $MODEL_NODES * $N_PROCS)"
      exit
   endif
endif

set n = 1
while ( $n <= $ENS_SIZE )

  if ( $RUN_IN_PBS == "yes" ) then  #  PBS queuing system

    echo "2i\"                                                                  >! advance.sed
    echo "#==================================================================\" >> advance.sed
    echo "#PBS -N init_mpas_${n}\"                                              >> advance.sed
    echo "#PBS -j oe\"                                                          >> advance.sed                                                             
    echo "#PBS -o logs/init_mpas_${n}.log\"                                     >> advance.sed
    echo "#PBS -A ${PROJ_NUMBER}\"                                              >> advance.sed
    echo "#PBS -q ${QUEUE}\"                                                    >> advance.sed
    echo "#PBS -l walltime=${TIME_INIT}\"                                       >> advance.sed
    echo "#PBS -l select=${MODEL_NODES}:ncpus=${N_CPUS}:mpiprocs=${N_PROCS}\"   >> advance.sed
    echo "#=================================================================="  >> advance.sed
    echo 's%${1}%'"${n}%g"                                                      >> advance.sed
    echo 's%${2}%'"${fn_param}%g"                                               >> advance.sed

    sed -f advance.sed ./init_mpas_grib.csh >! init_mpas_grib.pbs
    qsub init_mpas_grib.pbs
    ${REMOVE} init_mpas_grib.pbs advance.sed
    sleep 1

  else

    ./init_mpas_grib.csh $n $fn_param >! logs/init_mpas_${n}.log

  endif

  @ n++

end


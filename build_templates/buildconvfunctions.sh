#!/bin/bash

declare -a programs

source $DART/build_templates/buildpreprocess.sh

mpisrc="null_mpi"
windowsrc=""
m=""


#-------------------------
# print usage and exit
#-------------------------
function print_usage() {
  echo ""
  echo " Usage:   "
  echo "  buildconverter.sh               : build everything"
  echo "  buildconverter.sh clean         : clean the build" 
  echo "  buildconverter.sh help          : print help message"
  echo "   " 
  echo "  buildconverter.sh [program]     : optional argument " 
  echo "                                    [program] build a single program"
  echo "   " 
  exit
}


#-------------------------
# Remove programs, .o. .mod
#-------------------------
cleanup() {
\rm -f *.o *.mod Makefile

for p in ${programs[@]}; do  
  \rm -f $p
done
}

#-------------------------
# parse quicbuild.sh arguments:
#   nompi - compile without mpi
#   help - print usage and exit
#   program_name - compile single program
#   clean - remove programs
#-------------------------
function arguments() {

if [ $# -gt 2 ]; then
   print_usage
fi


case $1 in
  help)
    print_usage
    ;;

  clean)
    cleanup
    exit
    ;;
esac

single_prog=$1
}

#-------------------------
# Converter findsrc
#-------------------------
function findconvsrc() {

local core=$(find $DART/assimilation_code/modules -type d -name observations -prune -o -type f -name "*.f90" -print)
local conv=$(find $DART/observations/obs_converters/$CONVERTER -type f -name "*.f90" )
local convF90=$(find $DART/observations/obs_converters/$CONVERTER -type f -name "*.F90" )
local modelsrc="$DART/models/template/model_mod.f90"
local loc="$DART/assimilation_code/location/$LOCATION \
          $DART/assimilation_code/location/utilities/ \
          $DART/models/model_mod_tools/test_interpolate_$LOCATION.f90"
local misc="$DART/models/utilities/ \
            $DART/models/model_mod_tools/model_check_utilities_mod.f90 \
            $DART/observations/forward_operators/obs_def_mod.f90 \
            $DART//observations/forward_operators/obs_def_utilities_mod.f90 \
            $DART/assimilation_code/modules/observations/obs_kind_mod.f90 \
            $DART/assimilation_code/modules/observations/obs_sequence_mod.f90 \
            $DART/assimilation_code/modules/observations/forward_operator_mod.f90 \
            $DART/observations/obs_converters/utilities/obs_utilities_mod.f90"
local obserrsrc=$DART/observations/obs_converters/obs_error/$OBS_ERROR"_obs_err_mod.f90"

# remove null/mpi from list
mpi=$DART/assimilation_code/modules/utilities/mpi_utilities_mod.f90
nullmpi=$DART/assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
nullwin=$DART/assimilation_code/modules/utilities/null_win_mod.f90
craywin=$DART/assimilation_code/modules/utilities/cray_win_mod.f90
nocraywin=$DART/assimilation_code/modules/utilities/no_cray_win_mod.f90

if [ $mpisrc == "mpi" ]; then

   core=${core//$nullmpi/}
   core=${core//$nullwin/}
   if [ windowsrc == "craywin" ]; then
       core=${core//$nocraywin/}
   else #nocraywin
       core=${core//$craywin/}
   fi
else #nompi

   core=${core//$mpi/}
   core=${core//$nocraywin/}
   core=${core//$craywin/}
fi

convsrc="${core} ${conv} ${obserrsrc} ${modelsrc} ${misc} ${loc}"

# remove nuisance files
nuisance=(\
"$DART/assimilation_code/modules/assimilation/assim_tools_mod.pf.f90"
"$DART/assimilation_code/modules/assimilation/filter_mod.dopplerfold.f90"
"$DART/assimilation_code/modules/utilities/null_restart_pnetcdf_mod.f90"
"$DART/assimilation_code/modules/utilities/pnetcdf_utilities_mod.f90"
"$DART/assimilation_code/modules/utilities/restart_pnetcdf_mod.f90"
"$DART/models/bgrid_solo/fms_src/shared/time_manager/time_manager.f90"
)

for nus in  ${nuisance[@]}; do
 convsrc=${convsrc//$nus/}
done

# remove converter(s) {main.f90} from list

for p in ${programs[@]}; do
  prog=$DART/observations/obs_converters/$CONVERTER/$p.f90
  convsrc=${convsrc//$prog/}
done


}

#-------------------------
#-------------------------
function dartbuild() {

#look in $program directory for {main}.f90 
local program

if [ $1 == "obs_diag" ]; then
 program=$DART/assimilation_code/programs/obs_diag/$LOCATION
else
 program=$DART/assimilation_code/programs/$1
fi

 $DART/build_templates/mkmf -a $DART -x $m -p $1 \
     $EXTRA \
     $convsrc \
     $program \
     $DART/observations/obs_converters/$CONVERTER/$1.f90

}

#-------------------------
# buld converter
#-------------------------
buildconv() {

findconvsrc
i=1
for p in ${programs[@]}; do
  echo "Building " $p " build " $i " of " ${#programs[@]}
  dartbuild $p
  ((i++))
done


}


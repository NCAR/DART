#!/bin/bash

declare -a programs

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


# if the first argument is mpi or nompi
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
# Build and run preprocess
# Arguements: 
#  none
# Globals:
#  DART - root of DART
#-------------------------
function buildpreprocess() {

 # run preprocess if it is in the current directory
 if [ -f preprocess ]; then
   echo "already there"
   ./preprocess
   return
 fi

# link to preprocess if it is already built
if [ -f $DART/src/programs/preprocess/preprocess ]; then
   echo "not there but built"
   ln -s $DART/src/programs/preprocess/preprocess .
   ./preprocess 
   return
fi

 # build preproces
 cd $DART/src/programs/preprocess
 $DART/build_templates/mkmf -x -p $DART/src/programs/preprocess/preprocess \
      -a $DART $DART/src/programs/preprocess/path_names_preprocess
 cd -
 ln -s $DART/src/programs/preprocess/preprocess .
 ./preprocess
}


#-------------------------
# Converter findsrc
#-------------------------
function findconvsrc() {

local core=$(find $DART/src/core -type f -name "*.f90") 
#local modelsrc=$(find $DART/models/$MODEL/src -type d -name programs -prune -o -type f -name "*.f90" -print)
local modelsrc="$DART/models/template/model_mod.f90"
local conv="$DART/observations/obs_converters/$CONVERTER"
local loc="$DART/src/location/$LOCATION \
           $DART/src/model_mod_tools/test_interpolate_$LOCATION.f90"
local misc="$DART/src/location/utilities \
            $DART/models/utilities/default_model_mod.f90 \
            $DART/observations/forward_operators/obs_def_mod.f90 \
            $DART/observations/forward_operators/obs_def_utilities_mod.f90 \
            $DART/observations/obs_converters/utilities"
local obserrsrc=$DART/observations/obs_converters/obs_error/$OBS_ERROR"_obs_err_mod.f90"
local dew="$DART/observations/obs_converters/obs_error/dewpoint_obs_err_mod.f90 \
           $DART/observations/obs_converters/MADIS/meteor_mod.f90"


mpi=$DART/src/$mpisrc
window=$DART/src/$windowsrc  # TODO be carefull of this for null window

local conv=$(find $DART/observations/obs_converters/$CONVERTER -type f -name "*.f90" ) 

# remove converter {main.f90} from list
for p in ${programs[@]}; do
  prog=$DART/observations/obs_converters/$CONVERTER/$p.f90
  conv=${conv//$prog/}
done

convsrc="${core} ${conv} ${obserrsrc} ${dew} ${modelsrc} ${misc} ${loc} ${mpi} ${window}"

}

#-------------------------
#-------------------------
function dartbuild() {

#look in $program directory for {main}.f90 
local program

if [ $1 == "obs_diag" ]; then
 echo "Doing obs_diag" 
 program=$DART/src/programs/obs_diag/$LOCATION
else
 program=$DART/src/programs/$1
fi

echo $1.f90
 $DART/build_templates/mkmf -x $m -p $1 \
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


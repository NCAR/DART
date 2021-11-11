#!/bin/bash

declare -a programs
declare -a serial_programs
declare -a model_programs
declare -a model_serial_programs

#-------------------------
# print usage and exit
#-------------------------
function print_usage() {
  echo ""
  echo " Usage:   "
  echo "  quickbuild.sh                     : build everything"
  echo "  quickbuild.sh clean               : clean the build" 
  echo "  quickbuild.sh help                : print help message"
  echo "   " 
  echo "  quickbuild.sh [nompi] [program]   : optional arguments " 
  echo "                                      [nompi] build without mpi"
  echo "                                      [program] build a single program"
  echo "   " 
  echo "  Example 1. Build filter without mpi:"
  echo "           quickbuild.sh nompi filter"
  echo "   " 
  echo "  Example 2. Build perfect_model_obs with mpi"
  echo "           quickbuild.sh perfect_model_obs"
  echo "   " 
  exit
}

cleanup() {
\rm -f *.o *.mod Makefile
all_programs=("${programs[@]}" "${model_programs[@]}" "${serial_programs[@]}" "${model_serial_programs[@]}")

for p in ${all_programs[@]}; do 
  \rm -f $p
done
}

#-------------------------
# parse quickbuild.sh arguments:
#   nompi - compile without mpi
#   help - print usage and exit
#   program_name - compile single program
#   clean - remove programs
#-------------------------
function arguments() {

if [ $# -gt 2 ]; then
   print_usage
fi

# Default to build with mpi
mpisrc=mpi
windowsrc=no_cray_win
m="-w" # mkmf wrapper arg

# if the first argument is mpi or nompi
case $1 in
  help)
    print_usage
    ;;

  nompi)
    mpisrc="null_mpi"
    windowsrc=""
    m=""
    shift 1
    ;;

  clean)
    cleanup
    exit
    ;;
esac

single_prog=$1
}

#-------------------------
# Collect the source files needed to build DART
# Globals:
#  dartsrc - created buy this function
#  LOCATION - stores the location module
#  DART - expected in the enviroment
#  mpisrc  - quickbuild argument
#  windosrc - just no_cray win 
#-------------------------
function findsrc() {

local core=$(find $DART/src/core -type f -name "*.f90") 
local modelsrc=$(find $DART/models/$MODEL/src -type d -name programs -prune -o -type f -name "*.f90" -print)
local loc="$DART/src/location/$LOCATION \
           $DART/src/model_mod_tools/test_interpolate_$LOCATION.f90"
local misc="$DART/src/location/utilities \
            $DART/models/utilities/default_model_mod.f90 \
            $DART/observations/forward_operators/obs_def_mod.f90 \
            $DART/observations/forward_operators/obs_def_utilities_mod.f90"

mpi=$DART/src/$mpisrc
window=$DART/src/$windowsrc

dartsrc="${core} ${modelsrc} ${misc} ${loc} ${mpi} ${window}"
}

#-------------------------
# Build a program 
# Arguements: 
#  program name
# Globals:
#  DART - root of DART
#  dartsrc - source files
#  LOCATION - location mod (3D sphere, oned, etc.)
#  m - mpi mkmf wrapper
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

 $DART/build_templates/mkmf -x $m -p $1 \
     $dartsrc \
     $program
}

#-------------------------
# Build a model specific program
# looks in $DART/models/$MODEL/src/programs for {main}.f90 
# Arguements: 
#  program name
#  mpi wrapper arg for mkmf
# Globals:
#  DART - root of DART
#  dartsrc - source files
#-------------------------
function modelbuild() {
 $DART/build_templates/mkmf -x $m -p $1 $DART/models/$MODEL/src/programs/$1.f90 \
     $dartsrc
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
# Build executables
# Globals:
#  programs - array of DART programs
#  model_programs - array of model specific programs 
#  single_prog - name of a single program to build
#-------------------------
function buildit() {

if [ ! -z "$single_prog" ] ; then # build a single program
    if [[ " ${programs[*]} " =~ " ${single_prog} " ]]; then
       echo "building dart program " $single_prog
       findsrc
       dartbuild $single_prog
       exit
    elif [[ " ${serial_programs[*]} " =~ " ${single_prog} " ]]; then
       echo "building serial dart program " $single_prog
       mpisrc="null_mpi"
       windowsrc=""
       m=""
       findsrc
       dartbuild $single_prog
       exit
    elif [[ " ${model_programs[*]} " =~ " ${single_prog} " ]];then 
       echo "building model program" $single_prog
       findsrc
       modelbuild $single_prog
       exit
    elif [[ " ${model_serial_programs[*]} " =~ " ${single_prog} " ]];then 
       echo "building model program" $single_prog
       mpisrc="null_mpi"
       windowsrc=""
       m=""
       findsrc
       modelbuild $single_prog
       exit
    else
       echo "ERROR: unknown program" $single_prog
       exit 4
    fi  
fi

# if no single program argument, build everything
n=$((${#programs[@]}+${#model_programs[@]}+${#serial_programs[@]}+${#model_serial_programs[@]}))

i=1

# MPI programs 
findsrc
for p in ${programs[@]}; do
  echo "Building " $p " build " $i " of " $n
  dartbuild $p 
  ((i++))
done

for p in ${model_programs[@]}; do
  echo "Building " $p " build " $i " of " $n
  modelbuild $p 
  ((i++))
done

# clean before building serial programs
# if using mpi
[ $mpisrc == "mpi" ] && \rm -f *.o *.mod

mpisrc="null_mpi"
windowsrc=""
m=""

# Serial programs
findsrc

for p in ${serial_programs[@]}; do
  echo "Building " $p " build " $i " of " $n
  dartbuild $p  
  ((i++))
done

for p in ${model_serial_programs[@]}; do
  echo "Building " $p " build " $i " of " $n
  modelbuild $p  
  ((i++))
done

}

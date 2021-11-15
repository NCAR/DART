#!/bin/bash 

declare -a programs
declare -a serial_programs
declare -a model_programs
declare -a model_serial_programs

source $DART/build_templates/buildpreprocess.sh

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

local core=$(find $DART/assimilation_code/modules -type d -name observations -prune -o -type f -name "*.f90" -print)
local loc="$DART/assimilation_code/location/$LOCATION \
          $DART/assimilation_code/location/utilities/ \
          $DART/models/model_mod_tools/test_interpolate_$LOCATION.f90"
local modelsrc=$(find $DART/models/$MODEL -type d -name $EXCLUDE -prune -o -type f -name "*.f90" -print)
local misc="$DART/models/utilities/ \
            $DART/models/model_mod_tools/model_check_utilities_mod.f90 \
            $DART/observations/forward_operators/obs_def_mod.f90 \
            $DART//observations/forward_operators/obs_def_utilities_mod.f90 \
            $DART/assimilation_code/modules/observations/obs_kind_mod.f90 \
            $DART/assimilation_code/modules/observations/obs_sequence_mod.f90 \
            $DART/assimilation_code/modules/observations/forward_operator_mod.f90 \
            $DART/observations/obs_converters/utilities/obs_utilities_mod.f90"

# The quantity_mod.f90 files are in assimilation_code/modules/observations
# so adding individual files from assimilation_code/modules/observations

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

dartsrc="${core} ${modelsrc} ${loc} ${misc}"

# remove model specific programs from source list
all_modprogs=("${model_programs[@]}" "${model_serial_programs[@]}")
for modprog in ${all_modprogs[@]}; do
 mp=$DART/models/$MODEL/$modprog".f90"
 dartsrc=${dartsrc//$mp/}
done

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
 dartsrc=${dartsrc//$nus/}
done

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
 program=$DART/assimilation_code/programs/obs_diag/$LOCATION
else
 program=$DART/assimilation_code/programs/$1
fi

 $DART/build_templates/mkmf -x -a $DART $m -p $1 \
     $dartsrc \
     $EXTRA \
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
 $DART/build_templates/mkmf -x -a $DART $m -p $(basename $1) $DART/models/$MODEL/$1.f90 \
     $EXTRA \
     $dartsrc
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

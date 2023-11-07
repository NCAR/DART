#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#-------------------------
# Build functions for DART
#-------------------------
# Globals:
#  DART - path to DART directory
#       - expected in the enviroment
#  LOCATION - location module to use
#             set by quickbuild.sh
#  EXTRA - other source files that
#          are not in $MODEL directory
#          for example, files from another
#          model diretory.
#          set by quickbuild.sh
#  EXCLUDE - directories in model dir to exclude
#            from the collection of .f90 files
#            set by quickbuild.sh
#  MODEL - model being compiled, set by quickbuild.sh
#  dartsrc - created by findsrc
#  programs, serial_programs, model_programs, model_serial_programs
#          - arrays of program names to build
#            set by quickbuild.sh
#
#  dev_test - developer_test compilation. set by quickbuild.sh
#  TEST - used when compiling developer_tests. set by quickbuild.sh
#-------------------------

set -e
EXTRA=""
EXCLUDE=""
TEST=""
dev_test=0

declare -a programs
declare -a serial_programs
declare -a model_programs
declare -a model_serial_programs

source "$DART"/build_templates/buildpreprocess.sh

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
  echo "  quickbuild.sh [mpi/nompi/mpif08] [program]   : optional arguments " 
  echo "                                                 [mpi]     build with mpi (default)"
  echo "                                                 [nompi]   build without mpi"
  echo "                                                 [mpif08]  build with mpi using mpi_f08"
  echo "                                                 [program] build a single program"
  echo "   " 
  echo "  Example 1. Build filter without mpi:"
  echo "           quickbuild.sh nompi filter"
  echo "   " 
  echo "  Example 2. Build perfect_model_obs with mpi"
  echo "           quickbuild.sh perfect_model_obs"
  echo "   " 
  echo "  Example 3. Build perfect_model_obs with mpi using the mpi_f08 bindings"
  echo "           quickbuild.sh mpif08 perfect_model_obs"
  echo "   " 
  exit
}

#--------------------------
# Remove programs, .o. .mod
#--------------------------
cleanup() {

\rm -f -- *.o *.mod Makefile
\rm -f -- input.nml.*_default
all_programs=("${programs[@]}" "${model_programs[@]}" "${serial_programs[@]}" "${model_serial_programs[@]}")

for p in ${all_programs[@]}; do 
  \rm -f -- $(basename $p)
done

\rm -f -- preprocess
cleanpreprocess

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

# Default to build with mpi  (non f08 version)
mpisrc=mpi
windowsrc=no_cray_win
m="-w" # mkmf wrapper arg

# if the first argument is help, nompi, mpi, mpif08, clean
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

  mpi)
    mpisrc="mpi"
    windowsrc="no_cray_win"
    shift 1
    ;;  

  mpif08)
    mpisrc="mpif08"
    windowsrc="no_cray_winf08"
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
            $DART/assimilation_code/modules/observations/forward_operator_mod.f90"

# The quantity_mod.f90 files are in assimilation_code/modules/observations
# so adding individual files from assimilation_code/modules/observations

# remove null/mpi from list
local mpi="$DART"/assimilation_code/modules/utilities/mpi_utilities_mod.f90
local mpif08="$DART"/assimilation_code/modules/utilities/mpif08_utilities_mod.f90
local nullmpi="$DART"/assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
local nullwin="$DART"/assimilation_code/modules/utilities/null_win_mod.f90
local craywin="$DART"/assimilation_code/modules/utilities/cray_win_mod.f90
local nocraywin="$DART"/assimilation_code/modules/utilities/no_cray_win_mod.f90
local no_cray_winf08="$DART"/assimilation_code/modules/utilities/no_cray_winf08_mod.f90

if [ "$mpisrc" == "mpi" ]; then

   core=${core//$nullmpi/}
   core=${core//$nullwin/}
   core=${core//$mpif08/}
   core=${core//$no_cray_winf08/}
   if [ "$windowsrc" == "craywin" ]; then
       core=${core//$nocraywin/}
   else #nocraywin
       core=${core//$craywin/}
   fi

elif [ "$mpisrc" == "mpif08" ]; then

   core=${core//$nullmpi/}
   core=${core//$nullwin/}
   core=${core//$mpi/}
   core=${core//$craywin/}
   core=${core//$nocraywin/}

else  #nompi

   core=${core//$mpi/}
   core=${core//$mpif08/}
   core=${core//$nocraywin/}
   core=${core//$no_cray_winf08/}
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

for nus in  "${nuisance[@]}"; do
 dartsrc=${dartsrc//$nus/}
done

}

#-------------------------
# Build a program 
# Arguments: 
#  program name
#-------------------------
function dartbuild() {

local program

if [ $dev_test -eq 0 ]; then
  #look in $program directory for {main}.f90
  if [ $1 == "obs_diag" ]; then
    program=$DART/assimilation_code/programs/obs_diag/$LOCATION
  elif [ $1 == "streamflow_obs_diag" ]; then
    program=$DART/assimilation_code/programs/obs_diag/streamflow
  else
    program=$DART/assimilation_code/programs/$1
  fi
else
  # For developer tests {main}.f90 is in developer_tests
  program=$DART/developer_tests/$TEST/$1.f90
fi

 $DART/build_templates/mkmf -x -a $DART $m -p $1 \
     $dartsrc \
     $EXTRA \
     $program
}

#-------------------------
# Build a model specific program
# looks in $DART/models/$MODEL/src/programs for {main}.f90 
# Arguments:
#  program name
#-------------------------
function modelbuild() {
 $DART/build_templates/mkmf -x -a $DART $m -p $(basename $1) $DART/models/$MODEL/$1.f90 \
     $EXTRA \
     $dartsrc
}

#-------------------------
# Build executables
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

# when building all programs, remove input.nml.*_default files
\rm -f -- input.nml.*_default

}

#-------------------------
# Build any NetCDF files from .cdl files
#-------------------------
function cdl_to_netcdf {

shopt -s nullglob
for f in *.cdl ; do
  local outname=$(basename "$f" .cdl ).nc
  if [ ! -f "$outname" ]; then
    echo "building" "$outname" "from" "$f"
    ncgen -o "$outname" "$f"
  fi
done
shopt -u nullglob

}

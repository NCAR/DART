#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#-------------------------
# Build functions for observation converters
#-------------------------
# Globals:
#  DART - path to DART directory
#         expected in the enviroment
#  LOCATION - location module to use
#             set by quickbuild.sh
#  LIBRARIES - additional libraries to link to
#              for example libprepbufr.a
#              set by quickbuild.sh
#  EXTRA - other source files that
#          are not in $CONVERTER
#          for example, a specific model_mod
#          set by quickbuild.sh
#
#  convsrc - created by findconvsrc
#  programs - array of programs names to build
#             set by quickbuild.sh
#
# The GSI obs converter needs mpi
#  mpisrc="null_mpi"
#  windowsrc=""
#  m=""
#-------------------------
set -e
declare -a programs
source "$DART"/build_templates/buildpreprocess.sh

# Defaults
mpisrc="null_mpi"
windowsrc=""
m=""
LIBRARIES=""
EXTRA=""

#-------------------------
# print usage and exit
#-------------------------
function print_usage() {
  echo ""
  echo " Usage:   "
  echo "  quickbuild.sh               : build everything"
  echo "  quickbuild.sh clean         : clean the build"
  echo "  quickbuild.sh help          : print help message"
  echo "   " 
  echo "  quickbuild.sh [program]     : build a single program"
  echo "   " 
  exit
}


#--------------------------
# Remove programs, .o. .mod
#--------------------------
cleanup() {
\rm -f -- *.o *.mod Makefile input.nml*_default

for p in "${programs[@]}"; do
  \rm -f -- "$p"
done

\rm -f -- preprocess
cleanpreprocess

}

#-------------------------
# parse quickbuild.sh arguments:
#   help - print usage and exit
#   program_name - compile single program
#   clean - remove programs, .o .mod makefile
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
# Collect src needed to build converter
#-------------------------
function findconvsrc() {

local core=$(find "$DART"/assimilation_code/modules -type d -name observations -prune -o -type f -name "*.f90" -print)
local conv=$(find "$DART"/observations/obs_converters/"$CONVERTER" -type f -name "*.f90" )
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
else #nompi

   core=${core//$mpi/}
   core=${core//$mpif08/}
   core=${core//$nocraywin/}
   core=${core//$no_cray_winf08/}
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

for nus in  "${nuisance[@]}"; do
 convsrc=${convsrc//$nus/}
done

# remove converter(s) {main.f90} from list

for p in "${programs[@]}"; do
  prog="$DART"/observations/obs_converters/"$CONVERTER"/"$p".f90
  convsrc=${convsrc//$prog/}
done

}

#-------------------------
# Build a program 
# Arguments:
#  program name
#-------------------------
function dartbuild() {

#look in $program directory for {main}.f90 
local program
local mkmf_libs

if [ $1 == "obs_diag" ]; then
 program=$DART/assimilation_code/programs/obs_diag/$LOCATION
else
 program=$DART/assimilation_code/programs/$1
fi

# Additional libraries
if [ ! -z $LIBRARIES ]; then
  mkmf_libs="-l $LIBRARIES"
else
  mkmf_libs=""
fi
 
 $DART/build_templates/mkmf -a $DART -x $m $mkmf_libs -p $(basename $1)  \
     $EXTRA \
     $convsrc \
     $program \
     $DART/observations/obs_converters/$CONVERTER/$1.f90

}

#-------------------------
# build converter programs
#-------------------------
buildconv() {

findconvsrc

if [ ! -z "$single_prog" ] ; then # build a single program
    if [[ " ${programs[*]} " =~ " ${single_prog} " ]]; then
       echo "building dart program " $single_prog
       dartbuild $single_prog
       exit
    else
       echo "ERROR: unknown program " $single_prog
       exit 4
    fi
fi

# if no single program argument, build everything
i=1
for p in ${programs[@]}; do
  echo "Building " $p " build " $i " of " ${#programs[@]}
  dartbuild $p
  ((i++))
done

# when building all programs, remove input.nml.*_default files
\rm -f -- input.nml.*_default

}


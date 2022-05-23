#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-------------------------
# Build and run preprocess
# Arguements: 
#  none
# Globals:
#  DART - root of DART
#-------------------------
function buildpreprocess() {

 local pp_dir=$DART/assimilation_code/programs/preprocess

 # run preprocess if it is in the current directory
 if [ -f preprocess ]; then
   ./preprocess
   return
 fi

# link to preprocess if it is already built, run
if [ -f $pp_dir/preprocess ]; then
   ln -s $pp_dir/preprocess .
   ./preprocess 
   return
fi

 # build preprocess, link, run
 cd $pp_dir
 $DART/build_templates/mkmf -x -p $pp_dir/preprocess \
      -a $DART $pp_dir/path_names_preprocess
 cd -
 ln -s $pp_dir/preprocess .
 ./preprocess
}

#-------------------------
# clean up *.mod *.o for preprocess
#-------------------------
function cleanpreprocess() {

 local pp_dir=$DART/assimilation_code/programs/preprocess
 cd $pp_dir
 \rm -f -- *.o *.mod Makefile preprocess
 cd -

}


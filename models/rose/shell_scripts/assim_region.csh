#!/bin/csh

# set echo verbose

set PBS_O_WORKDIR = $1
set element = $2
set temp_dir = $3

# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

mkdir $temp_dir
cd $temp_dir

cp $PBS_O_WORKDIR/filter_assim_region__in$element filter_assim_region_in
cp ${PBS_O_WORKDIR}/assim_region .
cp ${PBS_O_WORKDIR}/filter_assim_obs_seq .

echo ls $temp_dir for element $element
echo junk > element$element
ls -lRt 

cp $PBS_O_WORKDIR/input.nml input.nml
cp $PBS_O_WORKDIR/rose.nml rose.nml

   ./assim_region > rose_out_temp
   echo after executing assim_region

echo element $element >> rose_reg_temp
ls -lt >> rose_reg_temp
cat rose_reg_temp >> $PBS_O_WORKDIR/rose_reg_temp$element
mv filter_assim_region_out $PBS_O_WORKDIR/filter_assim_region_out$element

cd $PBS_O_WORKDIR
rm -rf $temp_dir

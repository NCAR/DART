#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-- This is the real model advancement script called by advance_model.csh

set perfect_model = ${PERFECT}

set      process = $MP_CHILD
set   num_states = $1
set control_file = $2

set temp_dir = 'advance_temp'${process}

# Create a clean temporary directory and go there
\rm -rf  $temp_dir  || exit 1
mkdir -p $temp_dir  || exit 1
cd       $temp_dir  || exit 1

# Get the program and input.nml
ln  -s  ../input.nml        .  || exit 1

# Loop through each state
set state_copy = 1
set ensemble_member_line = `expr 3 \* $MP_CHILD + 1`
set      input_file_line = `expr $ensemble_member_line + 1`
set     output_file_line = `expr $ensemble_member_line + 2`

set ensemble_member = `head -$ensemble_member_line ../$control_file | tail -1`
set input_file      = `head -$input_file_line      ../$control_file | tail -1`
set output_file     = `head -$output_file_line     ../$control_file | tail -1`
   
#-- Get model.in from filter (input vector) for dynamo
ln -s ../$input_file ./temp_ic || exit 2
../dart_to_model

#-- If not the 1st step restart file should be available
if ( -f ../inrestart.dat.$ensemble_member ) then
   ln -s ../inrestart.dat.$ensemble_member ./inrestart.dat
endif 

#-- set model initialization data
ln -s ../inithi.dat .

if ( $perfect_model == 1 ) then

#-- Get the input vector for the day
   set day = `awk  'NF==2{print $2}' model.in|tail -1l`
   set phi = `../../perfect_model.sh $day`
   @ day--
   echo "0    $day" > model.out
   @ day++
   echo "0    $day" >> model.out
   set ct = 0
   while ( $ct < 1 )
      echo $phi >> model.out
      @ ct++
   end
   mv model.out model.in
#-- Model reads the input vector in model.in and advances
#-- to for some time
   ../dyn.out

#--store the model restart file
   cp outrestart.dat ../inrestart.dat.$ensemble_member

#-- prepare the input for model_to_dart
   echo "0    $day" > model.out
   awk 'NF==1{print $1}' model.in >> model.out

#-- Store observation for model_interpolate
   mv obsval.dat ../obsval.dat.$ensemble_member

else  #-- Assimilation steps

   ../dyn.out     #-- Run the dynamo model till next step --
   mv outrestart.dat outrestart.dat.bck
   cp obsval.dat obsval.dat.bck

   cp model.in model.in.bck #-- copy the input vector and then
   ../x_plus_delta          #-- change the vector to a small amount

#-- Run the dynamo with changed input vector but identical restart file
   ../dyn.out    

   ../dobsdx  #-- Calculate first derivative of observation, store in obsval.dat
   cp obsval.dat ../obsval.dat.$ensemble_member #-- for model_interpolate

   cp outrestart.dat.bck ../inrestart.dat.$ensemble_member #-- for next step

   cp model.in.bck model.in #-- restore the input vector as given by filter
   ../gau_param             #-- add little Gaussian noise to the input vector
endif
../model_to_dart            #-- Give it back to filter
mv temp_ud ../$output_file || exit 2

cd ..

\rm -rf $temp_dir

exit 0



#!/usr/bin/env bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# runs t2o multiple times to convert multiple CHAMP files and append the 
# resulting files to obs_seq.out (or the obs_out_file specified in 
# input.nml&text_to_obs_nml )

# grep gets the line, awk extracts the name
name=`grep "obs_out_file" < input.nml | awk -F\" '{for(i=2;i<=NF;i+=2)print $i}'`
echo $name
# remove the old file before starting the append-process
rm $name

# what is the first day in 2002 (year is hardcoded for now)
d1=335
# how many days do you want to append together?
# (note, if you need more than 2, you need to download them from
# http://sisko.colorado.edu/sutton/data/ver2.2/champ/density/2002/ascii/ )
nd=4

for (( i = 1 ; i <= $nd; i++ ))
do
    echo $i
    sed -i '.tmp' 's/text_input_file.*/text_input_file = "Density_3deg_02_'$[$d1+$i-1]'.ascii"/' input.nml
    ./text_to_obs > temp
done

#shows how many observations are there total
tail temp

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


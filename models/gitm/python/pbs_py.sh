#!/bin/sh
#
# DART $Id$
# from Alexey Morozov
#
#PBS -S /bin/sh
#PBS -N dart_test_tec_py
#PBS -l nodes=1:ppn=1,pmem=1000mb,walltime=30:00
#PBS -A ACCOUNT
#PBS -l qos=ACCOUNT
#PBS -q QUE
#PBS -M USERNAME@EXAMPLE.com
#PBS -m abe
#PBS -V

hd=~/HEAD/DART/models/gitm #directory where the "clean" (not temporary) gitm folder is located and where there is DART/observations folder. Only edit the python files (like plot_tec.py) in here, as these py file will overwrite the py files in $wd directory below
wd=/nobackup/USERNAME/dart_test_tec/gitm #directory where you want stuff plotted (it will look in ${wd}/work for *.nc files)
fn_gitm=${wd}/work/advance_temp_e21/UA #what truth-gitm run data files you want to plot (don't have to be in ${wd}) (as in 3DALL_t*.bin)
fn_gitm_type=3DALL #type of files provided above (only 3DALL is supported so far)
fn_gps=${hd}/../../observations/gnd_gps_vtec/work/gps021201g.002.txt #where the gps madrigal observation text file is (as in gps021201g.002.txt)
sim=False #was this gps file simulated (True) or real (False)? Affects what gets comparred to what in plot_tec.py and what gets plotted

cd $hd/python  || exit 1
cp $hd/python/*.py $wd/python  || exit 1 #in case I edited the py files 
cp $hd/python/*.txt $wd/python  || exit 1
cd $wd/python || exit 2

echo 'pbs_py.sh: past cp, starting pGITM'

cd $fn_gitm || exit 3
./pGITM #combines all block data files (*.b0001) into a one file defined over the whole globe (*.bin)
cd -

echo 'pbs_py.sh: past pGITM, starting python'

python plot_tec.py $fn_gitm/data/${fn_gitm_type}_t $fn_gps $sim || exit 3

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


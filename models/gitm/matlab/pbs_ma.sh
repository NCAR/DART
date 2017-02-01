#!/bin/sh
#
# DART $Id$
#
#PBS -S /bin/sh
#PBS -N dart_test_tec_ma
#PBS -l nodes=1:ppn=1,pmem=1000mb,walltime=30:00
#PBS -A ACCOUNT
#PBS -l qos=ACCOUNT
#PBS -q QUE
#PBS -M USERNAME@EXAMPLE.com
#PBS -m abe
#PBS -V

module load matlab

hd=~/HEAD/DART/models/gitm #directory where the "clean" (not temporary) gitm folder is located and where there is DART/observations folder. Only edit the python files (like plot_tec.py) in here, as these py file will overwrite the py files in $wd directory below
wd=/nobackup/NYXUSERNAME/dart_test_tec/gitm #directory where you want stuff plotted (it will look in ${wd}/work for *.nc files)
pbs_file_d=pbs_file.sh #name of the pbs file relative to $wd, which containes DART run information
truth_run_gitm=${wd}/work/advance_temp_e21/UA #what truth-gitm run data files you want to plot (don't have to be in ${wd}) (as in 3DALL_t*.b0001)
pbs_file_t=../../pbs_file.sh #name of the pbs file relative to $truth_run_gitm, which containes grid information (like nblockslon and nlons) and awt
middle_run_gitm=${wd}/work/advance_temp_e21/UA #what middle-gitm run data files you want to plot (don't have to be in ${wd}) (as in 3DALL_t*.b0001) #the middle run is different from the true one because of f107 it uses - I want to see what "middle of ensembe" would do if I just let it run without adjusting it. For example, if truth is 150 and ensemble is centered about 130, then I want to run three simulations: the truth (without assimilation, f=150), the middle (without assimilation, f=130) and the ensemble with assimilation (f is adjusted on the fly). right now middle_run_gitm is set to be the same as truth_run_dir, but that's only to check if all the files are in the right places - change this when you do the actual runs.
pbs_file_m=../../pbs_file.sh #name of the pbs file relative to $middle_run_gitm, which containes grid information (like nblockslon and nlons) and awt

cd $hd/matlab  || exit 1
cp $hd/matlab/*.m $wd/matlab  || exit 1 #in case I edited the .m files 
cd $wd/matlab || exit 2

echo 'pbs_ma.sh: past cp, starting pGITM'

cd $truth_run_gitm || exit 3
./pGITM  #combines all block data files (*.b0001) into a one file defined over the whole globe (*.bin)
cd -
cd $middle_run_gitm || exit 4
./pGITM  #combines all block data files (*.b0001) into a one file defined over the whole globe (*.bin)
cd -

echo 'pbs_ma.sh: past pGITM, starting matlab'

matlab -nodesktop -r "plot_champ('"$hd"/../..','"$truth_run_gitm"','"$pbs_file_t"','"$middle_run_gitm"','"$pbs_file_m"','"$wd/work"','"$pbs_file_d"')" || exit 5

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


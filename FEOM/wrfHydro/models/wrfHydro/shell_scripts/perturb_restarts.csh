#!/bin/csh

# Prupose: Given a WRF Hydro state and optionally a file of additional parameters to include
# in the assimilation state vector, create an initial ensemble. 

# Note: Perturbing the ensemble will likely dis-equilibrate the model and it therefore probably
# needs run forward in time before starting assimilation. See run_ensemble.csh for this. 

# Specifics:
# This script reads perturb.namelist for file and variable information. It also checks 
# the information against input.nml and hydro.namelist - so these two must be in accordance 
# with perturb.namelist before this script is run. 
#
# An assumption is that you have exactly two restart files matching both
# RESTART.*  &  HYDRO_RST.*
# in the directory 
# initialEnsembleDir 
# which is specified in perturb.namelist. 
# Optionally, a third could be included (but I havent worked out details of this yet):
# ASSIM_ONLY_RST.*
# These get converted to 
# restart.iEns.nc, restart.hydro.iEns.nc, & restart.noAssim.iEns.nc , respectively,
# where iEns represents the ith ensemble member. 
# *** These ensemble-specific files are always deleted if they exist when this 
# *** script is run.

#===============================================================================

# path to the perturb.namelist, should eventually just be "perturb.namelist"
set pertNml = 'perturb.namelist'

# get the number of ensembles.
set nEns = `grep -i 'nEns' $pertNml | tr -d ' '| egrep -v '^[!\]'`
set nEns = `echo $nEns | cut -d '=' -f2 | cut -d '!' -f1 | tr -d '! ,'`
# check the number of ensemble members in perturb.namelist against what's
# specified in input.nml.
## I decided this was a bad idea.
#set nEns2 = `grep -A 42 filter_nml input.nml | grep -i 'ens_size' | tr -d ' '| egrep -v '^[!\]'`
#set nEns2 = `echo ${nEns2} | cut -d '=' -f2 | cut -d '!' -f1 | tr -d '! ,'`
#if ( ${nEns} != ${nEns2} ) then
#    echo $nEns and $nEns2
#    echo "The number of ensembles specified in $pertNml and in input.nml DO NOT match."
#    exit 1
#endif 



# path to the geo fine grid / fulldom file. 
# This should match that in the hydro.namelist. 
set geoFineFile = `grep -i 'GEO_FINEGRID_FLNM' $pertNml | tr -d ' '| egrep -v '^[!\]'`
set geoFineFile = `echo $geoFineFile | cut -d '=' -f2 | cut -d '!' -f1 | tr -d "! ,\'"`
set geoFineFile2 = `grep -i 'GEO_FINEGRID_FLNM' hydro.namelist | tr -d ' "'| egrep -v '^[!\]'`
set geoFineFile2 = `echo $geoFineFile2 | cut -d '=' -f2 | cut -d '!' -f1 | tr -d "! ,\'"`

if ( ${geoFineFile} != ${geoFineFile2} ) then
    echo "The GEO_FINEGRID_FILENAME specified in $pertNml and in input.nml DO NOT match."
    exit 2
endif 

# get the path for the perturbed restart files. 
set pertPath = `grep -i 'initialEnsembleDir' $pertNml | tr -d ' '| egrep -v '^[!\]'`
set pertPath = `echo $pertPath | cut -d '=' -f2 | cut -d '!' -f1 | tr -d "! ,\'"`
## copy it over for posteritry
\cp -v --remove-destination $pertNml $pertPath/.

## the pert path needs to be an absolute path for now... so just automate 
set origDir = `pwd`
cd $pertPath
set fullPertPath = `pwd`
cd $origDir
ex perturb.namelist <<ex_end
g;initialEnsembleDir;s;= .*;= '$fullPertPath';
wq
ex_end
set pertPath = `grep -i 'initialEnsembleDir' $pertNml | tr -d ' '| egrep -v '^[!\]'`
set pertPath = `echo $pertPath | cut -d '=' -f2 | cut -d '!' -f1 | tr -d "! ,\'"`


## check that there is exactly one of each of LSM and Hydro restart files. 
set nLsmRestarts   = `ls -1   ${pertPath}/RESTART.* | wc -l`
if ( $nLsmRestarts != 1 ) then
    echo "There must be exactly one RESTART.* restart file in $pertPath"
    exit 3
endif
set nHydroRestarts = `ls -1 ${pertPath}/HYDRO_RST.* | wc -l`
if ( $nHydroRestarts != 1 ) then
    echo "There must be exactly one HYDRO_RST.* restart file in $pertPath"
    exit 4
endif
set nAssimOnlyRestarts = `ls -1 ${pertPath}/ASSIM_ONLY_RST.* | wc -l`
if ( $nAssimOnlyRestarts > 1 ) then
    echo "There must be no more than one ASSIM_ONLY_RST.* restart file in $pertPath"
    exit 5
endif

# clear any ensemble files (assuming ensemble sizes are < 1000)
\rm -f ${pertPath}/restart*0*.nc

echo $pertPath
# make fresh copes of the base restart files. 
foreach i (`seq 1 $nEns`)
    set ii=`printf "%04d" $i`
    \ls ${pertPath}/RESTART.* 
    \cp -v ${pertPath}/RESTART.* ${pertPath}/restart.${ii}.nc || exit 4
    \cp -v ${pertPath}/HYDRO_RST.* ${pertPath}/restart.hydro.${ii}.nc || exit 3
    if ( $nAssimOnlyRestarts > 0 ) then 
	\cp -v ${pertPath}/ASSIM_ONLY_RST.* ${pertPath}/restart.assimOnly.${ii}.nc || exit 4
    endif 
end

# perturb the restarts!
/home/jamesmcc/fortranTools/wrfHydroInitEns/smcIcs/smcPerturb.exe

echo '========================================================================================='
ls $pertPath
echo

exit 0


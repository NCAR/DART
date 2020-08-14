#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# ex: setup_experiment.csh forceCopyDartBuilds forceCopyDartScripts forceCopyModelParams forceCopyAll

#==============================================================================
# SSSS the commands so we can avoid problems with aliases, etc.

set   MOVE = '/bin/mv -v'
set   COPY = '/bin/cp -vp'
set   LINK = '/bin/ln -vs'
set REMOVE = '/bin/rm -rf'

#==============================================================================
## File structure initialization.

# where this script is run and where we look for and run the model
set CaseDirectory = `pwd`
echo 
echo 'CaseDirectory =' $CaseDirectory

## where we look for dart. Put the following in a ~/.wrf_hydro_tools file
set whtConfigFile=~/.wrf_hydro_tools   ## please dont commit changes
set whtFileMessage="\e[31mPlease create a $whtConfigFile file to contain wrf_hydro_dart=path/to/wrf_hydro_dart\e[0m"
if ( ! -e $whtConfigFile ) then
    echo "$whtFileMessage"
    exit 1
endif
set dartDir = `grep "wrf_hydro_dart=" $whtConfigFile | cut -d '=' -f2 | tr -d ' '` 
echo $dartDir
if ( $dartDir == '' ) then
    echo "$whtFileMessage"
    exit 1
endif

echo 'dartDir =' $dartDir

#==============================================================================
# copy this script, why not?

set setupFilterSource = ${dartDir}/models/wrf_hydro/shell_scripts/cheyenne/setup_experiment.csh 
set setupFilterTarget = setup_experiment.csh
if ( -e $setupFilterTarget ) then 
    if ( `diff -q $setupFilterSource $setupFilterTarget` !~ '' ) then 
	${COPY}   $setupFilterSource $setupFilterTarget
	echo "Updating this script: TRY AGAIN!"
	exit 2
    endif 
else 
    ${COPY} $setupFilterSource $setupFilterTarget || exit 2
endif 

#==============================================================================
## Name (logical) variables after input arguments, only modifies these 
## defaults if the proper variables are passed. Only need add a "var=0" line 
## for additional input args handled this way.
@ forceCopyDartBuilds = 0
@ forceCopyDartScripts = 0
@ forceCopyModelParams = 0
@ forceCopyAll = 0
@ n = 1
while ($n <= $#argv) 
    set $argv[$n] = 1
    @ n++
end

if ( $forceCopyAll ) then 
    @ forceCopyDartBuilds = 1
    @ forceCopyDartScripts = 1 
    @ forceCopyModelParams = 1
endif 
echo "forceCopyAll: $forceCopyAll"

#===============================================================================
# Clean up previous runs
# Destroy all previous output for the sake of clarity.
mv OUTPUT OUTPUT.DESTROYING
${REMOVE} OUTPUT.DESTROYING
mkdir OUTPUT

${REMOVE} member_????          
${REMOVE} dart_log.*             
${REMOVE} filter.log              
${REMOVE} filter_restart.*        
${REMOVE} restart*.nc             

#==============================================================================
# TJH: since input.nml is needed to determine active model components ... 
# this should come relatively early on.

echo
echo "DART Executables and run-time control"
foreach FILE ( filter perfect_model_obs input.nml )
    if ( -e ${FILE} && ! $forceCopyDartBuilds )  then
	echo "Using existing $FILE"
    else
	if(! -e ${dartDir}/models/wrf_hydro/work/${FILE} ) then 
	    echo "MISSING file: ${dartDir}/models/wrf_hydro/work/${FILE}"
	    exit 6
	endif
	${COPY} ${dartDir}/models/wrf_hydro/work/${FILE} . || exit 6
   endif
end

#==============================================================================
echo
echo "DART Scripts"
foreach FILE ( advance_ensemble.csh.template \
               run_filter.csh.template \
               submit_multiple_jobs_slurm.csh ) 
   if ( -e ${FILE} && ! $forceCopyDartScripts )  then
      echo "Using existing $FILE"
   else
      set FULLPATH = ${dartDir}/models/wrf_hydro/shell_scripts/cheyenne/${FILE}
	if (! -e ${FULLPATH}) then
	    echo "MISSING: ${FULLPATH}"
	    exit 7
	endif 
	${COPY} ${FULLPATH} . || exit 7
   endif
end

sed -e "s#EXPERIMENT_DIRECTORY#$CaseDirectory#" < \
       advance_ensemble.csh.template >! advance_ensemble.csh || exit 7
sed -e "s#EXPERIMENT_DIRECTORY#$CaseDirectory#" < \
       run_filter.csh.template       >! run_filter.csh || exit 7

chmod 755 advance_ensemble.csh run_filter.csh

#==============================================================================
# Check for required  files and dirs (excluding parameters)
foreach FILE ( wrf_hydro.exe namelist.hrldas hydro.namelist DOMAIN FORCING )
    if (! -e ${FILE} ) then
      echo "MISSING file: ${CaseDirectory}/${FILE}"
      exit 3
   endif
end

#===============================================================================
# Enforce some common sense namelist settings. 
# This might seem fussy but it will generally save alot of time in getting 
# everything set up correctly.

## require restart.nc
set MYSTRING = `grep RESTART_FILENAME_REQUESTED namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING = `echo $MYSTRING | cut -d'!' -f1 | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set restartFile = `echo $MYSTRING[$#MYSTRING]`
if ($restartFile !~ 'restart.nc' & $restartFile !~ './restart.nc' ) then
    echo "namelist.hrldas: The restart file (restart.nc) appears improperly specified in namelist.hrldas"
    exit 4
endif 

## require restart.hydro.nc
set MYSTRING = `grep RESTART_FILE hydro.namelist | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING = `echo $MYSTRING | cut -d'!' -f1 | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set hydroRestartFile = `echo $MYSTRING[$#MYSTRING]`
if ($hydroRestartFile !~ 'restart.hydro.nc' & $hydroRestartFile !~ './restart.hydro.nc' ) then
    echo "hydro.namelist: The restart file (restart.hydro.nc) appears improperly specified in hydro.namelist"
    exit 4
endif 

## require KHOUR = 1 
## I need to look into why this causes absolutely crazy behaviour where the 
## filter assimilates all obs at once and doesnt try to advance the mode.
set MYSTRING = `grep -i KHOUR namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING = `echo $MYSTRING | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set KHOUR = `echo $MYSTRING[$#MYSTRING]`
if ($KHOUR !~ '1') then
    echo "namelist.hrldas: KHOUR needs to be == 1 or you will hate life and filter will be crazy."
    exit 4
endif 

## require rst_dt = 1
#set MYSTRING = `grep -i rst_dt hydro.namelist | egrep -v '^[!]'`
set MYSTRING = `grep -i rst_dt hydro.namelist | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING = `echo $MYSTRING | cut -d'!' -f1 | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set rstDt = `echo $MYSTRING[$#MYSTRING]`
if ($rstDt !~ '60') then
    echo "hydro.namelist: rst_dt needs to be == 60."
    exit 4
endif 

# The following actually need to be tested from with a directory created at the current level
## might do something to make sure this doesnt/exist is unique... 

# WRFINPUT
set MYSTRING = `grep HRLDAS_SETUP_FILE namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING = `echo $MYSTRING | cut -d'!' -f1 | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set WRFINPUT = `echo $MYSTRING[$#MYSTRING]`

# FORCING 
set MYSTRING   = `grep INDIR namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING   = `echo $MYSTRING | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
set MYSTRING   = `echo $MYSTRING | sed -e 's#"# #g'`
set FORCINGDIR = `echo $MYSTRING[$#MYSTRING]`

# geo file for LSM
set MYSTRING   = `grep GEO_STATIC_FLNM namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING   = `echo $MYSTRING | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
set MYSTRING   = `echo $MYSTRING | sed -e 's#"# #g'`
set geoFile    = `echo $MYSTRING[$#MYSTRING]`

# geo file in hydro
set MYSTRING   = `grep GEO_STATIC_FLNM hydro.namelist | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING   = `echo $MYSTRING | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
set MYSTRING   = `echo $MYSTRING | sed -e 's#"# #g'`
set geoFileHydro = `echo $MYSTRING[$#MYSTRING]`

# Fulldom/finegrid file in hydro
set MYSTRING   = `grep GEO_FINEGRID_FLNM hydro.namelist | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING   = `echo $MYSTRING | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
set MYSTRING   = `echo $MYSTRING | sed -e 's#"# #g'`
set fineGridFile = `echo $MYSTRING[$#MYSTRING]`

## gw basin file, but only if it's set
set MYSTRING   = `grep -i GWBUCKPARM_file hydro.namelist | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING   = `echo $MYSTRING | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
set MYSTRING   = `echo $MYSTRING | sed -e 's#"# #g'`
set gwBasFile = `echo $MYSTRING[$#MYSTRING]`

set MYSTRING   = `grep GWBASESWCRT hydro.namelist | tr -d ' '| egrep -v '^[!\]'`
set MYSTRING   = `echo $MYSTRING | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
set MYSTRING   = `echo $MYSTRING | sed -e 's#"# #g'`
set gwBasSwitch = `echo $MYSTRING[$#MYSTRING]`

foreach file ($WRFINPUT $FORCINGDIR $geoFile $geoFileHydro \
              $fineGridFile $gwBasFile)
    echo "file: $file"
    if (! -e ${file}) then
 	if( $file == $gwBasFile & $gwBasSwitch == '0' ) continue
	echo $file 'Does NOT Exist as specified in a namelist'
	echo "Note that the path is relative to inside a directory within the calling level." 
	echo "You may need to make the symlink here (e.g: ln -s $file . )?"
	#cd ..
	#rm -rf tmpTestDir
	exit 4
    endif
end

#==============================================================================
# Parameters.
#   are gathered here for the model specfied in the namelist.hrldas
echo 
echo "PARAMS"
set paramsGathered = "PARAMS.gathered"
${REMOVE} $paramsGathered
\mkdir $paramsGathered

# Determine if this is a "noah" or a "noahmp" (lowercase)
# in order to copy the appropriate parameter files
# use the not-commented line specifying lsm_model_choice
set lsmModel = `grep -v '!' input.nml | grep -i lsm_model_choice`
set lsmModel = `echo $lsmModel | sed -e "s#[=,':]# #g"`
set lsmModel = `echo $lsmModel[$#lsmModel]`
set lsmModel = `echo $lsmModel | tr '[:upper:]' '[:lower:]'`

#noah params
## (the following three blocks should be written to use a function..)
if ( "$lsmModel" == "noah" ) then 
    set dir = "PARAMS.noah"
    if (! -e ${dir}) then 
	echo "MISSING directory in working dir: ${dir}/"
	exit 5
    endif 
    #foreach FILE (GENPARM.TBL LANDUSE.TBL SOILPARM.TBL URBPARM.TBL VEGPARM.TBL)
    foreach FILE (GENPARM.TBL SOILPARM.TBL URBPARM.TBL VEGPARM.TBL)
	if (! -e ${dir}/${FILE} )  then
	    echo "MISSING file: ${dir}/${FILE}"
	    exit 5
	endif
	echo "Copying ${dir}/${FILE} to $paramsGathered"
        ${COPY} ${dir}/${FILE} $paramsGathered || exit 5
    end
endif 

## noahMP params
if ( "$lsmModel" == "noahmp" ) then 
    set dir = "PARAMS.noahMP"
    if (! -e ${dir}) then 
	echo "MISSING directory: ${dir}"
	exit 5
    endif 
    foreach FILE (GENPARM.TBL MPTABLE.TBL SOILPARM.TBL URBPARM.TBL VEGPARM.TBL)
	if (! -e ${dir}/${FILE} )  then
	    echo "MISSING file: ${dir}/${FILE}"
	    exit 5
	endif
	echo "Copying ${dir}/${FILE} to $paramsGathered"
        ${COPY} ${dir}/${FILE} $paramsGathered || exit 5
    end
endif 
 
## hydro params
## (for now we'll assume the hydro component is always used. If it
##  is not being used, this is still probably harmless.)
set dir = "PARAMS.hydro"
if (! -e ${dir}) then 
    echo "MISSING directory: ${dir}"
    exit 5
endif 
foreach FILE ( CHANPARM.TBL HYDRO.TBL LAKEPARM.TBL )
    if (! -e ${dir}/${FILE} )  then
	echo "MISSING file: ${dir}/${FILE}"
	exit 5
    endif
    echo "Copying ${dir}/${FILE} to $paramsGathered"
    ${COPY} ${dir}/${FILE} $paramsGathered || exit 5
end 

#==============================================================================
# Initial ensemble

set ensembleDir = ${CaseDirectory}/initialEnsemble
echo
echo 'ensembleDir =' $ensembleDir

## if ensembleDir does not exist, create the directory but complain.
if (! -e ${ensembleDir} ) then 
    echo
    echo "CREATING $ensembleDir"
    \mkdir $ensembleDir  || exit 8
    echo "ENSEMBLE initial condition files have not been established (in the proper location)."
    echo "Please populate or link:  ${CaseDirectory}/${ensembleDir}"
    echo "mk_initial_ens.csh might help!"
    echo "run_initial_ens_fwd.csh might also help!"
    exit 8
endif 

ls -1 ${ensembleDir}/restart.[^ha]*         >! lsm_file_list.txt
ls -1 ${ensembleDir}/restart.hydro.*.nc     >! hydro_file_list.txt
ls -1 ${ensembleDir}/restart.assimOnly.*.nc >! param_file_list.txt

set   nLsmFiles = `cat   lsm_file_list.txt | wc -l`
set nHydroFiles = `cat hydro_file_list.txt | wc -l`
set nParamFiles = `cat param_file_list.txt | wc -l`

echo "Number of LSM   files: $nLsmFiles"
echo "Number of HYDRO files: $nHydroFiles"
echo "Number of PARAM files: $nParamFiles"

## Require one of lsm or hydro files
if ( $nLsmFiles < 1 & $nHydroFiles < 1 ) then
    echo "No initial restart files provided"
    echo "mk_initial_ens.csh might help!"
    echo "run_initial_ens_fwd.csh might also help!"
    exit 8
endif 

#==============================================================================
# DA inflation files
#
# A bit dodgy ... grab the variables that will be part of the DART state

set templateFile = `head -n 1 hydro_file_list.txt`
ncks -v qlink1 $templateFile temporary.nc
ncap2 -s 'qlink1(:)=1.2' temporary.nc input_priorinf_mean_d01.nc
ncap2 -s 'qlink1(:)=0.6' temporary.nc input_priorinf_sd_d01.nc





#==============================================================================

if ($nParamFiles > 0) then
    ${REMOVE} FORCING.perturbed || exit 9
    \mkdir    FORCING.perturbed || exit 9
    cd        FORCING.perturbed || exit 9
    foreach ii ( `seq $nParamFiles` )
      \mkdir ensemble.`printf "%04d" $ii`  || exit 9
    end
    cd .. || exit 9
endif

#==============================================================================
# Set the date in the namelist.hrldas to match that of the first restart file.

if ( $nLsmFiles > 0 ) then
    set restartFileTime = `ncdump -v Times $ensembleDir/restart.0001.nc | \
                           tail -n 2 | head -n 1 | cut -d'"' -f2`
else
    set restartFileTime = `ncdump -h $ensembleDir/restart.hydro.0001.nc | \
                           grep Restart_Time | cut -d'"' -f2`
endif
echo
echo "Setting namelist.hrldas restart time to match restarts in $ensembleDir"
echo
set restartFileYyyy = `echo $restartFileTime | cut -d- -f1`
set restartFileMm = `echo $restartFileTime | cut -d- -f2`
set restartFileDd = `echo $restartFileTime | cut -d- -f3 | cut -d_ -f1`
set restartFileHh = `echo $restartFileTime | cut -d_ -f2 | cut -d: -f1`

ex namelist.hrldas <<ex_end
g;START_YEAR;s;= .*;= $restartFileYyyy;
g;START_MONTH;s;= .*;= $restartFileMm;
g;START_DAY;s;= .*;= $restartFileDd;
g;START_HOUR;s;= .*;= $restartFileHh;
wq
ex_end

ex input.nml <<ex_end
/filter
g;ens_size ;s;= .*;= $nHydroFiles,;
g;num_output_state_members ;s;= .*;= $nHydroFiles,;
g;num_output_obs_members ;s;= .*;= $nHydroFiles,;
g;single_restart_file_in ;s;= .*;= .false.,;
g;single_restart_file_out ;s;= .*;= .false.,;
wq
ex_end

if (   $nLsmFiles < 1 ) then
   ${REMOVE} lsm_file_list.txt
endif
if ( $nHydroFiles < 1 ) then
   ${REMOVE} hydro_file_list.txt
endif
if ( $nParamFiles < 1 ) then
   ${REMOVE} param_file_list.txt
endif

#==============================================================================

exit 0


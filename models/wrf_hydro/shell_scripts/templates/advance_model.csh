#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.
# Can be used as-is with most low-order models and the bgrid model which
# can be advanced using the integrate_model executable.
# 
# Arguments are (created by 'filter' or 'perfect_model_obs' and include):
# 1) the process number of caller,
# 2) the number of ensemble members/state copies belonging to that process, and 
# 3) the name of the control_file for that process.
# 
# If this script finishes and the 'control_file' still exists, it is
# an ERROR CONDITION and means one or more of the ensemble members did
# not advance properly. Despite our best attempts to trap on this
# condition, some MPI installations simply hang, some properly terminate.
#
# This script loops over all the entries in the control_file to advance 
# any/all of the ensemble members.  The number of trips through the 
# loop is the second argument to this script. The control_file contains 
# the information about which ensemble members are to be advanced by THIS TASK.
# Sometimes it may be just one ensemble member, sometimes all of them.
# Read DART/doc/html/filter_async_modes.html and the mpi_intro.html
# for an overview.
#
# jlm An example layout of control and input files would be helpful. These links 
# dont provide much (quick) insight on that.
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance
#      and copies the necessary files into the temporary directory
# 2) copies/converts the DART state vector to something the model can ingest
#      (dart_to_wrf_hydro)
# 3) runs the model (wrf_hydro)
# 4) copies/converts the model output to input expected by DART
#      (wrf_hydro_to_dart)

set      process = $1
set   num_states = $2
set control_file = $3

set debug = `grep -i debug input.nml | cut -d'=' -f2 | cut -d',' -f 1 | tr -d ' '`
if ($debug > 0) echo "advance_model.csh: start"
#----------------------------------------------------------------------
# Block 1: copy necessary input files/executables/files common
#          to all model advances to a clean, temporary directory.
#          These will be used by ALL of the ensemble
#          members being advanced by this script.
#----------------------------------------------------------------------
if ($debug > 2) echo "advance_model.csh: block 1, setup files"

# Create a unique temporary working directory for this process's stuff
# The run-time directory for the entire experiment is called CENTRALDIR;
# we need to provide a safe haven for each TASK ... in 'temp_dir'.
set temp_dir = 'advance_temp'$process

set nodeDir  = `pwd`
if ( -e nodeDir.control ) then
    set nodeDir  = `head -1 nodeDir.control`
    if ( ! -e ${nodeDir}/OUTPUT ) mkdir -p ${nodeDir}/OUTPUT
    set temp_dir = ${nodeDir}/${temp_dir}
    set thisHost = `hostname`
    echo $thisHost >> nodeDir.control
endif

## this is used to remove relative paths 1 dir up.
## what happens if larger relative paths are used?
cd ..
set loginUp1Dir = `pwd`
cd -    
set loginDir = `pwd`

# Create a clean temporary directory and go there
\rm -rf   $temp_dir  || exit 1
\mkdir -p $temp_dir  || exit 1
cd        $temp_dir  || exit 1

# Get the DART input.nml and the NOAH namelist
\cp ${loginDir}/namelist.hrldas .
\cp ${loginDir}/hydro.namelist .
\cp ${loginDir}/input.nml .

# if running on node disks then have to modify the paths in the namelist files. 
##
#if ($?loginUp1Dir) then
    foreach ff ( namelist.hrldas hydro.namelist )
        if ($debug > 2) echo $ff
        foreach i ( `egrep '[.][.][/]' $ff | tr -d ' ' | egrep -v '^(!)' | cut -d'=' -f2` )
            if ($debug > 2) echo $i
            set whLine = `grep -n $i $ff | cut -d':' -f1`
            sed -i "${whLine}s&../&${loginUp1Dir}/&" $ff
        end
    end
#endif 

## Get the DOMAIN
ln -s ${loginDir}/DOMAIN .
ln -s ${loginDir}/FORCING .

    
# Get parameters.
foreach FILE ( ${loginDir}/PARAMS.gathered/* ) 
    if ($debug > 2) echo $FILE
#   \cp -v ${loginDir}/$FILE . || exit 2 ## if these are to be changed, that's handled below.
    \ln -sf $FILE . || exit 2  
end

## JLM TODO: figure if lsm_model_active and if hydro_model_active
set lsm_model_active1 = `grep -v '!' input.nml | grep -i lsm_model_choice | cut -d= -f2 | tr -cd '[[:alnum:]]._-' | wc -m`
set lsm_model_active = 0
if ($lsm_model_active > 0) set lsm_model_active = 1

set hydro_model_active1 = `grep -v '!' input.nml | grep -i hydro_model_active | cut -d= -f2 | tr -cd '[[:alnum:]]._-'`
set hydro_model_active = 0
if ($hydro_model_active1 == ".true.") set hydro_model_active = 1

# From the namelist determine if the noAssim restarts are needed. 
# The line could be commented out (default is blank in model_mod.f90) or set to ''.
set assimOnly_active1 = `grep -v '!' input.nml | grep -i assimOnly_netcdf_filename | wc -l`
set assimOnly_active2 = `grep -v '!' input.nml | grep -i assimOnly_state_variables | \
                         cut -d= -f2 | tr -cd '[[:alnum:]]._-' | wc -m`
set assimOnly_active = 0
if ($assimOnly_active1 > 0 & $assimOnly_active2 > 0) set assimOnly_active = 1
# assimOnly variables dont necessarily mean perturbed forcing?
# for now the do mean perturbed forcing...
# eventually may want to list the variables in the file and set flags based on these. 

#-------------------------------------------------------------------------------
# Loop through each state / ensemble member
set state_copy = 1
# These reference line number in the control file.
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)
   
    set ensemble_member = `\head -$ensemble_member_line ${loginDir}/$control_file | \tail -1`
    set input_file      = `\head -$input_file_line      ${loginDir}/$control_file | \tail -1`
    set output_file     = `\head -$output_file_line     ${loginDir}/$control_file | \tail -1`
    set instance        = `printf "%04d" $ensemble_member`

    #-------------------------------------------------------------------
    # Block 2: dart_to_wrf_hydro
    #          * remove scraps from previous advances
    #          * copy/link ensemble-member-specific files
    #          * convey the advance-to-time to the model
    #          * dart_to_wrf_hydro: convert the DART state vector to model format
    #-------------------------------------------------------------------

    if ($debug > 2) echo "advance_model.csh: block 2, converting ensemble member $instance"

    # clean up from last advance
    # some of these must be copied at some point?? for diagnostics?
    \rm -f  restart.nc  restart.hydro.nc  dart_restart  restart.assimOnly.nc
    \rm -f  wrf_hydro_advance_information.txt
    \rm -f  HYDRO_RST.*  RESTART.* 
    # if perturbed forcings are used, there will be *LDASIN_DOMAIN* here. see below.
    \rm -f  *.LDASOUT_DOMAIN*  *LDASIN_DOMAIN*
    \rm -f  *.LSMOUT_DOMAIN*  *.RTOUT_DOMAIN*  *.CHRTOUT*  *.CHANOBS*  frxst_pts_out.txt
    \rm -f  qstrmvol*  diag_hydro.*  stderr.txt stdout.txt  GW_*.txt  *.GW_DOMAIN*

    # need the wrf_hydro restart files for the output of dart_to_wrf_hydro
    if ($lsm_model_active > 0) \
        \ln -sv ${loginDir}/restart.$instance.nc  restart.nc   || exit 2
    if ($hydro_model_active > 0) \
        \ln -sv ${loginDir}/restart.hydro.$instance.nc  restart.hydro.nc   || exit 2

    if ( $assimOnly_active ) \
	\ln -sv ${loginDir}/restart.assimOnly.$instance.nc  restart.assimOnly.nc   || exit 2

    if ( $assimOnly_active ) then 
	## these guys depend on their initial values in call to apply_assimOnly
	## e.g. need to know the original geo_finegrid file to copy and modify it.
	\cp ${loginDir}/namelist.hrldas .
	\cp ${loginDir}/hydro.namelist .
	#\cp ${temp_dir}/namelist.hrldas .
	#\cp ${temp_dir}/hydro.namelist .
        # Since parameter files are potentially altered in assimOnly,
	# refresh and dont use symlinks.
	foreach FILE ( ${loginDir}/PARAMS.gathered/* ) 
	    echo $FILE
	    \cp -v --remove-destination $FILE . || exit 2 ## these might change so dont link
	end

    endif

    # the input file is the name of the dart_restart.instance? 
    # There must be reasons for being cryptic.
    \ln -sv ${loginDir}/$input_file           dart_restart || exit 2

    # push the assimilation back to the model
    # This modifies
    # restart.nc -> ${loginDir}/restart.$instance.nc
    # restart.hydro.nc -> ${loginDir}/restart.hydro.$instance.nc
    # restart.assimOnly.nc -> ${loginDir}/restart.assimOnly.$instance.nc
    echo `pwd`
    ${loginDir}/dart_to_wrf_hydro                          || exit 2

    if ( ! -e wrf_hydro_advance_information.txt ) then
	echo "ERROR: dart_to_noah failed for member $ensemble_member"
	echo "ERROR: dart_to_noah failed for member $ensemble_member"
	echo "ERROR: dart_to_noah failed for member $ensemble_member"
        exit 1
    endif

    # This next two parts are based on using one-hour forcing files
    # since the minimum time to advance the model seems to be 1 hour.
    # (kday, khour, but no kminute, for example)
    # dart_to_wrf_hydro provides the setting for namelist.hrldas:khour
    # we need to put that value in the local copy of namelist.hrldas

    set numadvancestr = `\grep -i khour wrf_hydro_advance_information.txt`
    set numadvancestr = `echo $numadvancestr | sed -e "s#[=,']# #g"`
    set numadvancestr = `echo $numadvancestr | sed -e 's#"# #g'`
    set numadvances   = `echo $numadvancestr[$#numadvancestr]`

# seems most efficient to only write restarts after the desired advance. 
# todo?

    ## ALSO have to keep the start time in hrldas abreast of the advancing.
    ## else the forcing data seems to have no effect.
    if ($lsm_model_active > 0) then
        set restartFileTime = `ncdump -v Times restart.nc | tail -2 | head -1 | cut -d'"' -f2`
    else
        set restartFileTime = `ncdump -h restart.hydro.nc | grep Restart_Time | cut -d'"' -f2`
    endif
    
    set restartFileYyyy = `echo $restartFileTime | cut -d- -f1`
    set restartFileMm = `echo $restartFileTime | cut -d- -f2`
    set restartFileDd = `echo $restartFileTime | cut -d- -f3 | cut -d_ -f1`
    set restartFileHh = `echo $restartFileTime | cut -d_ -f2 | cut -d: -f1`

ex namelist.hrldas <<ex_end
g;KHOUR;s;=.*;= $numadvances;
g;START_YEAR;s;=.*;= $restartFileYyyy;
g;START_MONTH;s;=.*;= $restartFileMm;
g;START_DAY;s;=.*;= $restartFileDd;
g;START_HOUR;s;=.*;= $restartFileHh;
wq
ex_end

    echo '******************************************************************************'
    echo ' Ensemble member:' $instance
    egrep 'START_(Y|M|D|H)' namelist.hrldas | grep -v !
    echo '******************************************************************************'

    # The forcing has to be for the NEXT "FORCING_TIMESTEP", apparently.
    # FORCING_TIMESTEP is defined in namelist.input At this point, dart_to_wrf_hydro
    # has assumptions that the forcing_timestep is one hour.

    # grep -n identifies the (line number): in the outupt, this becomes skipNlines
    set numfilestring = `\grep -ni nfiles wrf_hydro_advance_information.txt`
    set numfilestring = `echo $numfilestring | sed -e "s#[=,':]# #g"`
    set numfilestring = `echo $numfilestring | sed -e 's#"# #g'`
    set numfiles      = `echo $numfilestring[$#numfilestring]`
    set skipNlines    = `echo $numfilestring[1]`

    #-------------------------------------------------------------------
    # Block 2.5: AssimOnly Assimilation State Variables.
    if ( $assimOnly_active ) then 
      cp ${loginDir}/apply_assimOnly.csh .
      ./apply_assimOnly.csh ${loginDir}
      set failure = $?
      if ( $failure ) exit 22
    endif
    
    #-------------------------------------------------------------------
    # Block 3: advance the model
    #          In this case, we are saving the run-time messages to
    #          a LOCAL file, which makes debugging easier.
    #          integrate_model is hardcoded to expect input in temp_ic 
    #          and it creates temp_ud as output. 
    #          Your model will likely be different.
    #-------------------------------------------------------------------
    if ($debug > 2) echo "advance_model.csh: block 3, advance the model"
    if ($debug <= 2) echo "advance the model"
    
    set logTime=`date +%Y-%m-%d_%H.%M.%S`
    #mpirun -np 2 ${loginDir}/wrf_hydro.exe |& tee ${logTime}.stdout
    ## The next line is for shared nodes on cheyenne. It DOES use mpirun instead of mpiexec_mpt
    #mpirun `hostname` -np 2 ${loginDir}/wrf_hydro.exe >& ${logTime}.stdout 
    mpiexec_mpt -np 2 ${loginDir}/wrf_hydro.exe >& ${logTime}.stdout

    @ numadvancesNum = $numadvances

    if ($lsm_model_active) then
        @ lsm_status = `\ls -1 RESTART*DOMAIN* | wc -l`
        if ( $lsm_status < 1 )  then
            echo "ERROR: wrf_hydro died without RESTART files."
            \ls -l
            exit 23
        endif 
	if ( $lsm_status > $numadvancesNum ) then 
	    \ls -l RESTART*DOMAIN*
	    echo "WARNING: wrf_hydro created the above RESTART files. only expected # $numadvances" 
	endif 
    endif

    if ($hydro_model_active) then
        @ hydro_status = `\ls -1 HYDRO_RST* | wc -l`
        if ( $hydro_model_active > 0 && $hydro_status < 1 )  then
            echo "ERROR: wrf_hydro died without HYDRO_RST files."
            \ls -l
            exit 23
        endif
        if ( $hydro_status > $numadvances ) then 
	    \ls -l HYDRO_RST*
	    echo "WARNING: wrf_hydro created the above HYDRO_RST files. Only expected # $numadvances" 
	endif 
    endif


    #-------------------------------------------------------------------
    # Block 4: wrf_hydro_to_dart and managing output files. 
    #          rename files to reflect the ensemble member ID
    #-------------------------------------------------------------------
    if ($debug > 2) echo "advance_model.csh: block 4, manage/rename output files"
    # Do this before setting up the next run as there is an unwatned hydro restart file.

    # Determine model integration period of interest (which may contain multiple indiv
    # model integrations) and create a directory in central dir appropriately stamped 
    # to catch the output during this period.. 
    # The timestamps of the first and last LDASOUT files give us this even though we 
    # dont want the last LDASOUT file (it's for one hour beyond the desired integration period
    # thanks to HRLDAS).
    if ($lsm_model_active) then
        set integStart = `\ls -1 *LDASOUT_DOMAIN* | \head -n1 | \cut -d. -f1`
        set integEnd = `\ls -1 *LDASOUT_DOMAIN* | \tail -n1 | \cut -d. -f1`
    else
        set integStart = `\ls -1 *CHRTOUT_DOMAIN* | \head -n1 | \cut -d. -f1`
        set integEnd = `\ls -1 *CHRTOUT_DOMAIN* | \tail -n1 | \cut -d. -f1`
    endif

    set integStartYMD = `echo $integStart | cut -c1-8`
    set integStartH = `echo $integStart | cut -c9-10`
    set integStartM = `echo $integStart | cut -c11-12`
    set timeStep = `grep -i OUTPUT_TIMESTEP namelist.hrldas | cut -d'=' -f2 | tr -d ' '`
    set integStart = `date -d "$integStartYMD ${integStartH}:${integStartM} ${timeStep} seconds ago" +%Y%m%d%H%M`
        
    set integDir =  OUTPUT/model_integration.${integStart}-${integEnd}.$instance
    set integDirCurrent = ${nodeDir}/${integDir}
    \mkdir $integDirCurrent

    # Fixed HRLDAS to do restarts at the end of the loop, after time advance,  
    # with/after wrf_hydro. So I dont have to clean up a bunch of files.

    # Move the output files (*not* restarts)
    # Arezoo: Have removed these since it is not generated at this verison, and cause the script to fail ....
#    foreach outFile ( GW_*.txt frxst_pts_out.txt qstrmvolrt_accum.txt )
#	\mv $outFile ${integDirCurrent}/.
#    end 
    # these have their own timestamps but tag them with ensId/instance. 
    foreach outFile ( *.LDASOUT_DOMAIN* *.LSMOUT_DOMAIN* *.RTOUT_DOMAIN* *.CHRTOUT* *.CHANOBS* )
	\mv $outFile ${integDirCurrent}/${outFile}.${instance}.nc
    end 

   ## These are the parameters which were used in the model advance. 
   if ( $assimOnly_active ) cp restart.assimOnly.nc $integDirCurrent/.

    # Set the new/latest restart for ingest to dart
    if ($lsm_model_active) then
        set RESTARTlsm = `\ls -1  RESTART* | \tail -1`
        \ln -sf $RESTARTlsm    restart.nc        || exit 4
    endif
    
    if ($hydro_model_active) then
        set RESTARThydro = `\ls -1  HYDRO_RST* | \tail -1`
        \ln -sf $RESTARThydro  restart.hydro.nc  || exit 4
    endif

    ${loginDir}/wrf_hydro_to_dart              || exit 4

    \mv -v  dart_ics  ${loginDir}/$output_file          || exit 5
    # this breaks the restart.nc and restart.hydro.nc symlinks
    # but they are reset in the next loop
    if ($lsm_model_active) then
        \mv -v  ${RESTARTlsm}    ${integDirCurrent}/${RESTARTlsm}.$instance.nc   || exit 5
        \ln -sfv ${integDirCurrent}/${RESTARTlsm}.$instance.nc   ${loginDir}/restart.$instance.nc
    endif
    
    if ($hydro_model_active) then
        \mv -v  ${RESTARThydro}  ${integDirCurrent}/${RESTARThydro}.$instance.nc || exit 5
        \ln -sfv ${integDirCurrent}/${RESTARThydro}.$instance.nc ${loginDir}/restart.hydro.$instance.nc
    endif
    
    # the linking (vs. cp ing) in these last two lines implies that the model-created 
    # restart files are not sacred, they will be overwritten by dart_to_wrf_hydro

    ## increment
    @ state_copy++
    @ ensemble_member_line = $ensemble_member_line + 3
    @ input_file_line = $input_file_line + 3
    @ output_file_line = $output_file_line + 3

echo "advance_model.csh end loop over member for this process"

end  ## loop over ensemble members for this process.

# Change back to original directory and get rid of temporary directory.
# If all goes well, there should be no need to keep this directory.
# If you are debugging, you may want to keep this directory. 
cd $loginDir
\rm -rf $temp_dir

# MANDATORY - Remove the control_file to signal completion. If it still
# exists in CENTRALDIR after all the ensemble members have been advanced,
# it means one or more of the advances failed and is an ERROR CONDITION.

\rm -rf $control_file

## i think this code chunk is obselete.....
## if parameters in the geoFineFile were edited clean up that file. 
if ($?geoFineFileSrc) then
    echo 'fofofofofofofofofofofofofofofofofofofofofofofofofofofofof'
    \cp $geoFineFileOrig $geoFineFile
    \rm -f $geoFineFileOrig
endif

if ($debug > 0) echo "advance_model.csh: finish"

exit 0


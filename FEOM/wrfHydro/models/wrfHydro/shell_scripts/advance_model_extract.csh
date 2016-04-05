    if ( $assimOnly_active ) then 

	# Get a list of the variables in the restart.assimOnly file
	set assimOnlyVars = `ncks -m restart.assimOnly.nc | grep -E ': type' | cut -f 1 -d':'`

	#-----------------------------------------------------------
	# GEO_FINEGRID_FLNM parameters
	# if any vars which use this file, the original file needs to be copied 
	# and edits performed on the original file *location/name*. Note that 
	# the DOMAIN directory is a symlink, not the files in it. 
	# All of this gets cleaned up after the advances. 
	set geoFineFileParams = '(OVROUGHRTFAC|RETDEPRTFAC)'
	if ($instance == 0001) then 
	    if (`echo $assimOnlyVars | egrep -ci $geoFineFileParams`) then 
		# Copy and deal with the symlink that is/should be $geoFineFile
		set geoFineFile = `grep -v '!' hydro.namelist | grep -i GEO_FINEGRID_FLNM`
		set geoFineFile = `echo $geoFineFile | cut -d'"' -f2`
		set geoFineFileOrig = ${geoFineFile}.ORIG
		\cp $geoFineFile $geoFineFileOrig
	    endif
	endif 

	if (`echo $assimOnlyVars | egrep -ci $geoFineFileParams`) then 
	    rm  -f $geoFineFile
	    set operationString = ''
	    ## for now I'm only handling scalars. NCO commands should handle larger fields. 
	    ## but I might need to query the dimensions and/or deal with basin masks.

	    if (`echo $assimOnlyVars | egrep -ci OVROUGHRTFAC`) then 
		set ovRough = `ncks -H -v OVROUGHRTFAC ../restart.assimOnly.${instance}.nc \
				| head -1 | cut -d= -f2`
		set operationString = `echo $operationString "OVROUGHRTFAC=OVROUGHRTFAC*${ovRough}"`
	    endif 
	    
	    if (`echo $assimOnlyVars | egrep -ci RETDEPRTFAC`) then 
		set retDep = `ncks -H -v RETDEPRTFAC ../restart.assimOnly.${instance}.nc \
				| head -1 | cut -d= -f2`
		set operationString = `echo $operationString "RETDEPRTFAC=RETDEPRTFAC*${retDep}"`
	    endif 

	    ncap2 -s $operationString $geoFineFileOrig $geoFineFile

	endif 

	#-----------------------------------------------------------
	# Perturbed Forcing (just precip multiplier currently)
	# It may be desirable to keep the perturbed forcings (for diagnostics or smoother?).
	# Here you can exercise that choice. If kept, perturbed forcings will be stored to
	# ../FORCING.perturbed/yyyymmddhh.iEns.LDASIN_DOMAIN* where the final character 
	# is the same as in ../FORCING/yyyymmddhh.LDASIN_DOMAIN*. 
	# If cycling or re-running is performed for previous timesteps after each analysis, 
	# this may become more complex, with the forcing at a given time actually depending 
	# on the time of the analysis. (That is the forcing is changing as the analysis moves
	# into the future). But Im' not going to worry about that yet. 
	#
	# Create the forcings in the current directory and only move if they are to be kept.  
	# The FORGING.perturbed/* filenames should have date and the instance/nnnn should be specified.
	# then symlinked to remove this.
	## Alter the namelist.hrldas to point INDIR='./'
	if (`echo $assimOnlyVars | grep -ci precipMult`) then 

	    ## The name of the forcing file is conveniently supplied in 
	    ## wrfHydro_advance_information.txt. Though ( fix ) the times in the 2 top lines
	    ## seem wrong.
	
	    @ ifile = 1
	    while ($ifile <= $numfiles)
		@ linenum = $skipNlines + $ifile
		set FNAME = `\head -$linenum wrfHydro_advance_information.txt | tail -1`
		set FDATE = `echo $FNAME | sed -e "s#[.,']# #g"`
		set FDATE = `echo $FDATE[1]`

		set FFILE = `\ls ../FORCING/$FDATE.LDASIN_DOMAIN*`
		set FFILElocal = `echo $FFILE | \cut -d'/' -f3`
		
		echo $FFILE
		echo $FFILElocal
		    
		## get the precip multiplier out of the restart.assimOnly.nc file.
		## this will only work for scalar precip! 
		set precipMult = `ncks -H -v precipMult ../restart.assimOnly.${instance}.nc \
				| head -1 | cut -d= -f2`
				    
		## multiply and put the restulting file in the right place
		ncap2 -s "RAINRATE=RAINRATE*${precipMult}" \
			$FFILE ../FORCING.perturbed/ensemble.${instance}/$FFILElocal
			    
		## change the location of the input in the namelist.hrldas
		## maybe this could be done elsewhere, outside loop, but it's lightweight.
ex namelist.hrldas <<ex_end
g;INDIR;s;= .*;= "../FORCING.perturbed/ensemble.${instance}/";
wq
ex_end

		## need to update the timestamp in the restart.
		## the timestamp should match that of the forcing file??
		## jlm - fixme

		@ ifile = $ifile + 1
	    end
	endif # precip multiplier
	

    endif 

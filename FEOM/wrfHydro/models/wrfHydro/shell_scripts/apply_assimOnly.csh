#!/bin/csh

# Purpose: 
#  Apply forcing adjustments, and parameters adjustments (etc) supplied in 
#  restart.assimOnly.nc restart files to the model run in current directory. 

# Assumptions: 
#  You are in a directory with 
#  1) one restart.assimOnly.nc
#  2) all the other information necessary to advance the model (incl. parameter and 
#     geolocation files) are properly specified in the namelist.hrldas and 
#     hydro.namelist which also reside in this dir. A summary of 
#     file information used:
#     namelist.hrldas
#         KHOUR: the number of forcing timesteps to which forcing adjustments are applied. 
#         INDIR: where the forcing files are found
#     hydro.namelist
#         GEO_FINEGRID_FLNM: where certain model parameters are set

# Notes: 
#  1) Previously, when the model was to be advanced for long periods, I would 
#     consider the specified precipMult as a mean and I would generate a timeseries
#     about this mean with some variance. It's not clear that's appropriate. So it has
#     been abandoned. It seems like it could simply be added in at a given timestep
#     by introducing zero biased noise (rather than generating a full timeseries 
#     upfront).
#  2) Following on 1, eliminating the addtional perturbation eliminates any need to 
#     keep perturbed timeseries because they can be recreated from restart.assimOnly.nc
#     and the forcing. If the additional perturbation is desired, then the numbers/noise
#     should be just written to a separate file. 

#===============================================================================

set loginDir=`pwd`
if ( $?1 ) set loginDir=$1
echo 'loginDir: '$loginDir

## Setup
set restartTime = `ncdump -h restart.hydro.nc | grep Restart_ | cut -d'=' -f2 | tr -d ' ";'`
set startYyyy = `echo $restartTime | cut -d- -f1`
set startMm   = `echo $restartTime | cut -d- -f2`
set startDd   = `echo $restartTime | cut -d- -f3 | cut -d_ -f1`
set startHh   = `echo $restartTime | cut -d_ -f2 | cut -d: -f1`

## these need to be gathered from input.nml not restart.assimOnly
## that will provide greater flexibility with using subsets of the restart.assimOnly
##BAD: set assimOnlyVars = `ncks -m restart.assimOnly.nc | grep ': type' | cut -f1 -d':'`
## this identifies the line where the second name list starts
set whStartNml2 = `cat input.nml | tr -d ' ' | egrep -n '^[&]' | head -2 | tail -1 | cut -d':' -f1`
## this identifies the line where the assimOnly_state_variables start
set whStartAssimOnly = `cat input.nml | tr -d ' ' | egrep -n '^assimOnly_state_variables' | cut -d':' -f1`
if ( whStartAssimOnly == '' ) then
    exit 0
else 
    @ lengthAssimOnlyVars = $whStartNml2 - $whStartAssimOnly
    ## i'd like to break this up but handling of new line is annoying in csh (maybe elsewhere too), i'd have to add a "tr ' ' '\n'" 
    set assimOnlyVars = `tail input.nml -n+$whStartAssimOnly | head -$lengthAssimOnlyVars | tr -d "' " | tr -d '"' | egrep -v '^!' | egrep ',' | cut -d'=' -f2 | cut -d',' -f1`
endif 
#===============================================================================
# Forcings
# Edit the path to the forcing dir in the namelist.hrldas to use perturbed forcing
# in the current dir

set forcingVars = ( precipMult otherVarsHere )
set forcingGrep = `echo $forcingVars | tr ' ' '|'`
set forcingGrep = `echo "($forcingGrep)"`
set forcingActive =  `echo $assimOnlyVars | egrep $forcingGrep | wc -l`

# Original forcing location . Note that there no required name for this.
¡±# We need this even if forcing is not active because the path is one dir deeper
set inDirIn = `grep INDIR namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
set inDirIn = `echo $inDirIn | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
set inDirIn = `echo $inDirIn | sed -e 's#"# #g'`
set inDirIn = `echo $inDirIn[$#inDirIn]`

echo 'inDirIn: '$inDirIn

if ( $forcingActive ) then 

    set forcDir = FORCING.assimOnly
    ## break it down to build it back up (make sure it's clean)
    ## cant delete it in this routine since the model is called after exit.
    \rm -rf $forcDir
    \mkdir $forcDir

#make the namelist aware of the new forcings
ex namelist.hrldas <<ex_end
g;INDIR ;s;= .*;= "${forcDir}/";
wq
ex_end


    set forcType = `grep FORC_TYP namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
    set forcType = `echo $forcType | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
    set forcType = `echo $forcType | sed -e 's#"# #g'`
    set forcType = `echo $forcType[$#forcType]`
    if ( $forcType != 1 & $forcType != 6 ) then 
        echo "FORC_TYP = $forcType in namelist.hrldas NOT SUPPORTED CURRENTLY! (apply_assimOnly.csh)"
        exit 8
    endif

    # How many forcing times need adjusted?
    set khourIn = `grep KHOUR namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
    set khourIn = `echo $khourIn | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
    set khourIn = `echo $khourIn | sed -e 's#"# #g'`
    set khourIn = `echo $khourIn[$#khourIn]`

    set endYyyy = `date -ud "UTC $startYyyy-$startMm-$startDd $startHh hour + $khourIn hours" +%Y`

    ## This is not easy to do/optimize when there are lots of files. This is my best shot.
    ## It is relatively fast. 
    set yearSeq = `seq -s, $startYyyy $endYyyy | tr ',' '|'`
    ## the annchor (^) supposedly gives a speed up
    ## we always need these 
    set forcFilesSuper = \
        `find $inDirIn/ -regextype posix-extended \
            -regex "^$inDirIn.*/(${yearSeq}).*LDASIN.*" | sort`
    set forcFilesPosn = `echo $forcFilesSuper | tr ' ' '\n' | \
                          grep -n "${startYyyy}${startMm}${startDd}${startHh}" | cut -d':' -f1`
    set forcFilesIn = `echo $forcFilesSuper | tr ' ' '\n' | tail -n+$forcFilesPosn | head -$khourIn`
    \cp $forcFilesIn $forcDir/.
    set forcFiles = `ls $forcDir/*`

    if ( $forcType == 6 ) then 
        set forcFilesSuper6 = \
            `find $inDirIn/ -regextype posix-extended \
                -regex "^$inDirIn.*/(${yearSeq}).*PRECIP_FORCING.*" | sort`    
        set forcFilesPosn6 = `echo ${forcFilesSuper6} | tr ' ' '\n' | \
                                grep -n "${startYyyy}${startMm}${startDd}${startHh}" | cut -d':' -f1`
        set forcFilesIn6 = `echo ${forcFilesSuper6} | tr ' ' '\n' | tail -n+${forcFilesPosn6} | head -$khourIn`
        \cp $forcFilesIn6 $forcDir/.
        set forcFiles = `ls $forcDir/*PRECIP_FORCING*`
    endif            

    #-----------------------------------------------------------
    # precipMult
    # See if it is among the assimOnlyVars
    if ( `echo $assimOnlyVars | grep precipMult | wc -l` ) then 
        set precipMultDims = `ncks -m restart.assimOnly.nc | grep 'precipMult dimension' | wc -l`

        if ( $precipMultDims == 1 ) then
            set thePrecipMult = `ncdump -v precipMult restart.assimOnly.nc | tail -n2 | head -1 | \
                                    cut -d'=' -f2 | tr -d ' ;'`
            foreach iForc ( $forcFiles )
                if ( $forcType == 1 ) ncap2 -O -s "RAINRATE=RAINRATE*${thePrecipMult}" $iForc $iForc
                if ( $forcType == 6 ) ncap2 -O -s "precip_rate=precip_rate*${thePrecipMult}" $iForc $iForc
            end
        endif 

        if ( $precipMultDims == 2 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # precipMult

    #-----------------------------------------------------------
    # more vars... to come.

else  ## forcingActive else not
    ## we are one dir level deeper than we should be for the ensembles
    cd ../   
    if ( ! -e $inDirIn ) ln -sf $inDirIn .
    cd -
endif

#===============================================================================
# Parameters
set paramVars = ( OVROUGHRTFAC \
                  RETDEPRTFAC  \
                  gwCoeff      \
                  gwExpon      \
                  ksatMult     \
                  slope        \
                  refDk        \
                  refKdt       \
                  maxSmcMult   \
                  wltSmc       \
                  refSmc       \
                  satDw        \
                  rsMult       \
                  ch2opMult    \
                  czil         \
                  bbMult       \
                  satPsiMult   \
                  mann         \
)

set paramGrep = `echo $paramVars | tr ' ' '|'`
set paramGrep = `echo "($paramGrep)"`
set paramActive =  `echo $assimOnlyVars | egrep $paramGrep | wc -l`

if ( $paramActive ) then 

    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
    # GEO_FINEGRID: params which live in this file need set there 
    # for each ensemble member individually. copy to local and change the 
    # location of the file in use.
    set geoFineVars = ( OVROUGHRTFAC RETDEPRTFAC )
    set geoFineGrep = `echo $geoFineVars | tr ' ' '|'`
    set geoFineGrep = `echo "($geoFineGrep)"`
    set geoFineActive =  `echo $assimOnlyVars | egrep $geoFineGrep | wc -l`

    if ($geoFineActive) then 

        # original file
        set geoFineFileIn = `grep GEO_FINEGRID_FLNM hydro.namelist | tr -d ' '| egrep -v '^[!\]'`
        set geoFineFileIn = `echo $geoFineFileIn | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
        set geoFineFileIn = `echo $geoFineFileIn | sed -e 's#"# #g'`
        set geoFineFileIn = `echo $geoFineFileIn[$#geoFineFileIn]`
        # ensemble member specific file
        set geoFineFile   = `echo $geoFineFileIn | tr '/' ' '`
        set geoFineFile   = `echo $geoFineFile[$#geoFineFile]`
        ## copy the file  locally 
        \cp $geoFineFileIn $geoFineFile
        ## change location
ex hydro.namelist <<ex_end
g;GEO_FINEGRID_FLNM;s;=.*;= "$geoFineFile";
wq
ex_end

        # OVROUGHRTFAC -------
        # See if it is among the assimOnlyVars
        if ( `echo $assimOnlyVars | grep OVROUGHRTFAC | wc -l` ) then 

            set ovRoughDims = `ncks -m restart.assimOnly.nc | grep 'OVROUGHRTFAC dimension' | wc -l`

            if ( $ovRoughDims == 1 ) then
                set theOvRough = `ncdump -v OVROUGHRTFAC restart.assimOnly.nc | tail -n2 | head -1 | \
                                cut -d'=' -f2 | tr -d ' ;'`
                ncap2 -O -s "OVROUGHRTFAC=OVROUGHRTFAC*0+${theOvRough}f" $geoFineFile $geoFineFile
            endif 
                
            if ( $ovRoughDims == 2 ) then
                echo 'Not yet configured'
                exit 3
            endif 
        endif  # OVROUGHRTFAC
   

        # RETDEPRTFAC --------
        # See if it is among the assimOnlyVars
        if ( `echo $assimOnlyVars | grep RETDEPRTFAC | wc -l` ) then 

            set retDepDims = `ncks -m restart.assimOnly.nc | grep 'RETDEPRTFAC dimension' | wc -l`

            if ( $retDepDims == 1 ) then
                set theRetDep = `ncdump -v RETDEPRTFAC restart.assimOnly.nc | tail -n2 | head -1 | \
                                cut -d'=' -f2 | tr -d ' ;'`
                ncap2 -O -s "RETDEPRTFAC=RETDEPRTFAC*0+${theRetDep}f" $geoFineFile $geoFineFile
            endif 
                
            if ( $retDepDims == 2 ) then
                echo 'Not yet configured'
                exit 3
            endif 
        endif  # RETDEPRTFAC

    endif # geoFineActive
 


    # PARAMETER FILES - *.TBL
    # If running filter, these WERE ALREADY copied locally by advance_model.csh if the 
    # assimOnly is active. OTherwise, it dosent seem to matter as adjustments are to a 
    # tmp file which is *moved* to the original.
    # NOTE: the files are format sensitive.

    #--------------------------------------------    
    #--------------------------------------------    
    # CHANPARM.TBL 
    
    # mann
    if ( `echo $assimOnlyVars | grep mann | wc -l` ) then 
        set mannDims = `ncks -m restart.assimOnly.nc | grep 'mann dimension' | wc -l`
        if ( $mannDims == 1 ) then
            set theMann = `ncdump -v mann restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            set theMann = `echo $theMann | cut -c1-6`
            \rm -rf CHANPARM.TBL.NEW
            touch   CHANPARM.TBL.NEW
            set nlines = `cat CHANPARM.TBL | wc -l`
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3)$' | wc -l` ) then
                    sed -n ${ll}p CHANPARM.TBL >> CHANPARM.TBL.NEW || exit 7
                    continue
                endif
                set origValue = `sed -n ${ll}p CHANPARM.TBL | cut -d',' -f5`
                set origValue = `printf '%.16f' $origValue`
                set newValue  = 0`echo "$origValue * $theMann" | bc | cut -c1-3`
                sed -n ${ll}p CHANPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$5="  "newValue;print}' \
                    >> CHANPARM.TBL.NEW 
            end

            mv CHANPARM.TBL.NEW   CHANPARM.TBL
        endif 
        if ( $mannDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  


    #--------------------------------------------    
    #--------------------------------------------    
    ## GWBUCKPARM.TBL

    ## gwCoeff ---------------
    if ( `echo $assimOnlyVars | grep gwCoeff | wc -l` ) then 
        set gwCoeffDims = `ncks -m restart.assimOnly.nc | grep 'gwCoeff dimension' | wc -l`
        if ( $gwCoeffDims == 1 ) then
            set theGwCoeff = `ncdump -v gwCoeff restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            set theGwCoeff = `echo $theGwCoeff | cut -c1-6`
            cat GWBUCKPARM.TBL | \
                awk -v theGwCoeff=$theGwCoeff 'BEGIN{FS=",";OFS=","};{$2=theGwCoeff;print}' \
                > GWBUCKPARM.TBL.NEW
            mv GWBUCKPARM.TBL.NEW GWBUCKPARM.TBL
        endif 
        if ( $gwCoeffDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 

    ## gwExpon ---------------
    if ( `echo $assimOnlyVars | grep gwExpon | wc -l` ) then 
        set gwExponDims = `ncks -m restart.assimOnly.nc | grep 'gwExpon dimension' | wc -l`
        if ( $gwExponDims == 1 ) then
            set theGwExpon = `ncdump -v gwExpon restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            set theGwExpon = `echo $theGwExpon | cut -c1-6`
            \rm -rf GWBUCKPARM.TBL.NEW
            touch GWBUCKPARM.TBL.NEW
            sed -n 1p GWBUCKPARM.TBL > GWBUCKPARM.TBL.NEW
            sed -n 2p GWBUCKPARM.TBL | \
                awk -v theGwExpon=$theGwExpon 'BEGIN{FS=",";OFS=","};{$3=theGwExpon;print}' \
                >> GWBUCKPARM.TBL.NEW
            mv GWBUCKPARM.TBL.NEW GWBUCKPARM.TBL
        endif 
        if ( $gwExponDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 

    #--------------------------------------------    
    #--------------------------------------------    
    ## SOILPARM.TBL and HYDRO.TBL
    ## SOILPARM code must be above GENPARM
    ## Oh joy, some of these parameters these are in two places (SOILPARM and HYDRO!

    ## ksatMult --------------
    ## an example of multiplying a column
    if ( `echo $assimOnlyVars | grep ksatMult | wc -l` ) then 

        set ksatMultDims = `ncks -m restart.assimOnly.nc | grep 'ksatMult dimension' | wc -l`
        if ( $ksatMultDims == 1 ) then
            set theKsatMult = `ncdump -v ksatMult restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## SOILPARM.TBL
            set nlines = `cat SOILPARM.TBL | wc -l`
            \rm -rf SOILPARM.TBL.NEW
            touch SOILPARM.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3|23|24|25|45)$' | wc -l` ) then
                    sed -n ${ll}p SOILPARM.TBL >> SOILPARM.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f8`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theKsatMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p SOILPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$8="  "newValue;print}' \
                    >> SOILPARM.TBL.NEW 
            end

            mv SOILPARM.TBL.NEW   SOILPARM.TBL

            ## HYDRO.TBL
            set nlines = `cat HYDRO.TBL | wc -l`
            \rm -rf HYDRO.TBL.NEW
            touch HYDRO.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( $ll <= 31 ) then
                    sed -n ${ll}p HYDRO.TBL >> HYDRO.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p HYDRO.TBL | cut -d',' -f1`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theKsatMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p HYDRO.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$1=newValue;print}' \
                    >> HYDRO.TBL.NEW 
            end

            mv HYDRO.TBL.NEW   HYDRO.TBL

        endif ## 1D

        if ( $ksatMultDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 

    ## maxSmcMult ------------
    if ( `echo $assimOnlyVars | grep maxSmcMult | wc -l` ) then 

        set maxSmcMultDims = `ncks -m restart.assimOnly.nc | grep 'maxSmcMult dimension' | wc -l`
        if ( $maxSmcMultDims == 1 ) then
            set theMaxSmcMult = `ncdump -v maxSmcMult restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## SOILPARM.TBL
            set nlines = `cat SOILPARM.TBL | wc -l`
            \rm -rf SOILPARM.TBL.NEW
            touch SOILPARM.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3|23|24|25|45)$' | wc -l` ) then
                    sed -n ${ll}p SOILPARM.TBL >> SOILPARM.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f5`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theMaxSmcMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p SOILPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$5="  "newValue;print}' \
                    >> SOILPARM.TBL.NEW 
            end

            mv SOILPARM.TBL.NEW   SOILPARM.TBL

            ## HYDRO.TBL
            set nlines = `cat HYDRO.TBL | wc -l`
            \rm -rf HYDRO.TBL.NEW
            touch HYDRO.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( $ll <= 31 ) then
                    sed -n ${ll}p HYDRO.TBL >> HYDRO.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p HYDRO.TBL | cut -d',' -f2`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theMaxSmcMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p HYDRO.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$2=newValue;print}' \
                    >> HYDRO.TBL.NEW 
            end

            mv HYDRO.TBL.NEW   HYDRO.TBL

        endif ## 1D

        if ( $maxSmcMultDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 

echo foo12 >&2
    ## bbMult ------------
    if ( `echo $assimOnlyVars | grep bbMult | wc -l` ) then 

        set bbMultDims = `ncks -m restart.assimOnly.nc | grep 'bbMult dimension' | wc -l`
        if ( $bbMultDims == 1 ) then
            set theBbMult = `ncdump -v bbMult restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## SOILPARM.TBL
            set nlines = `cat SOILPARM.TBL | wc -l`
            \rm -rf SOILPARM.TBL.NEW
            touch SOILPARM.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3|23|24|25|45)$' | wc -l` ) then
                    sed -n ${ll}p SOILPARM.TBL >> SOILPARM.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f2`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theBbMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p SOILPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$2="  "newValue;print}' \
                    >> SOILPARM.TBL.NEW 
            end

            mv SOILPARM.TBL.NEW   SOILPARM.TBL

        endif ## 1D

        if ( $bbMultDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 

    endif  # bbMult

    ## satPsiMult ------------
    if ( `echo $assimOnlyVars | grep satPsiMult | wc -l` ) then 

        set satPsiMultDims = `ncks -m restart.assimOnly.nc | grep 'satPsiMult dimension' | wc -l`
        if ( $satPsiMultDims == 1 ) then
            set theSatPsiMult = `ncdump -v satPsiMult restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## SOILPARM.TBL
            set nlines = `cat SOILPARM.TBL | wc -l`
            \rm -rf SOILPARM.TBL.NEW
            touch SOILPARM.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3|23|24|25|45)$' | wc -l` ) then
                    sed -n ${ll}p SOILPARM.TBL >> SOILPARM.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f7`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theSatPsiMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p SOILPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$7="  "newValue;print}' \
                    >> SOILPARM.TBL.NEW 
            end

            mv SOILPARM.TBL.NEW   SOILPARM.TBL

        endif ## 1D

        if ( $satPsiMultDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 

    endif  # satPsiMult


    
    ## refSmc ------------
    ## This must be below maxsmc, satdk and bb
    ## if set to -1, set using the Campbell '74 equations suggested in the NoahMP code
    if ( `echo $assimOnlyVars | grep refSmc | wc -l` ) then 

        set refSmcDims = `ncks -m restart.assimOnly.nc | grep 'refSmc dimension' | wc -l`
        if ( $refSmcDims == 1 ) then
            set theRefSmc = `ncdump -v refSmc restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## SOILPARM.TBL
            set nlines = `cat SOILPARM.TBL | wc -l`
            \rm -rf SOILPARM.TBL.NEW
            touch SOILPARM.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3|23|24|25|45)$' | wc -l` ) then
                    sed -n ${ll}p SOILPARM.TBL >> SOILPARM.TBL.NEW || exit 7
                    continue
                endif

                ## if calculating from other table values, need maxsmc, satdk, and bb from table
                if ( `echo "$theRefSmc == -1" | bc` ) then                
                    if ( `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f12 | tr -d "'"` == "WATER" ) then
                        set newValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f6`
                    else
                        set tblMaxSmc = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f5 | sed -e 's/[eE]+*/\\*10\\^/'`
                        set tblSatDk  = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f8 | sed -e 's/[eE]+*/\\*10\\^/'`
                        set tblBb     = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f2 | sed -e 's/[eE]+*/\\*10\\^/'`
                        ## REFSMC1=MAXSMC*(5.79E-9/SATDK)**(1/(2*BB+3)) 5.79E-9 m/s= 0.5 mm
                        set theBase = `echo "scale=16; ((5.79*10^-9)/($tblSatDk))"  | bc`
                        set theExp  = `echo "scale=16; 1./(2.*$tblBb+3.)" | bc`
                        set thePow  = `echo "scale=16; e($theExp*l($theBase))" | bc -l`
                        set newValue  = `echo "scale=16; $tblMaxSmc*$thePow" | bc`
                    endif 
                else 
                    set origValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f6`
                    set origValue = `printf '%.16f' $origValue`
                    set newValue = `echo "$origValue * $theRefSmc" | bc`
                endif
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p SOILPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$6="  "newValue;print}' \
                    >> SOILPARM.TBL.NEW 
            end

            mv SOILPARM.TBL.NEW   SOILPARM.TBL
        endif ## 1D

        if ( $refSmcDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 

    endif  # refSmc


    ## wltSmc ------------
    ## This must be below maxsmc, satdk and bb
    ## if set to -1, set using the Campbell '74 equations suggested in the NoahMP code
    if ( `echo $assimOnlyVars | grep wltSmc | wc -l` ) then 

        set wltSmcDims = `ncks -m restart.assimOnly.nc | grep 'wltSmc dimension' | wc -l`
        if ( $wltSmcDims == 1 ) then
            set theWltSmc = `ncdump -v wltSmc restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## SOILPARM.TBL
            set nlines = `cat SOILPARM.TBL | wc -l`
            \rm -rf SOILPARM.TBL.NEW
            touch SOILPARM.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3|23|24|25|45)$' | wc -l` ) then
                    sed -n ${ll}p SOILPARM.TBL >> SOILPARM.TBL.NEW || exit 7
                    continue
                endif

                ## if calculating from other table values, need maxsmc, satdk, and bb from table
                if ( `echo "$theWltSmc == -1" | bc` ) then                
                    if ( `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f12 | tr -d "'"` == "WATER" ) then
                        set newValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f10`
                    else
                        set tblMaxSmc = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f5 | sed -e 's/[eE]+*/\\*10\\^/'`
                        set tblSatPsi  = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f7 | sed -e 's/[eE]+*/\\*10\\^/'`
                        set tblBb     = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f2 | sed -e 's/[eE]+*/\\*10\\^/'`
                        ## WLTSMC1=MAXSMC*(200./SATPSI)**(-1./BB) 
                        set theBase = `echo "scale=16; 200./$tblSatPsi"  | bc`
                        set theExp  = `echo "scale=16; -1/$tblBb" | bc`
                        set thePow  = `echo "scale=16; e($theExp*l($theBase))" | bc -l`
                        set newValue = `echo "scale=16; $tblMaxSmc*$thePow" | bc`
                    endif 
                else 
                    set origValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f10`
                    set origValue = `printf '%.16f' $origValue`
                    set newValue = `echo "$origValue * $theWltSmc" | bc`
                endif
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p SOILPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$10="  "newValue;print}' \
                    >> SOILPARM.TBL.NEW 
            end

            mv SOILPARM.TBL.NEW   SOILPARM.TBL
        endif ## 1D

        if ( $wltSmcDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 

    endif  # wltSmc


    ## satDw ------------
    ## This must be below maxsmc, satdk and bb
    ## if set to -1, set using the Campbell '74 equations suggested in the NoahMP code
    if ( `echo $assimOnlyVars | grep satDw | wc -l` ) then 

        set satDwDims = `ncks -m restart.assimOnly.nc | grep 'satDw dimension' | wc -l`
        if ( $satDwDims == 1 ) then
            set theSatDw = `ncdump -v satDw restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## SOILPARM.TBL
            set nlines = `cat SOILPARM.TBL | wc -l`
            \rm -rf SOILPARM.TBL.NEW
            touch SOILPARM.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3|23|24|25|45)$' | wc -l` ) then
                    sed -n ${ll}p SOILPARM.TBL >> SOILPARM.TBL.NEW || exit 7
                    continue
                endif

                ## if calculating from other table values, need maxsmc, satdk, and bb from table
                if ( `echo "$theSatDw == -1" | bc` ) then                
                    if ( `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f12 | tr -d "'"` == "WATER" ) then
                        set newValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f9`
                    else
                        set tblMaxSmc = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f5 | sed -e 's/[eE]+*/\\*10\\^/'`
                        set tblSatDk  = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f8 | sed -e 's/[eE]+*/\\*10\\^/'`
                        set tblSatPsi  = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f7 | sed -e 's/[eE]+*/\\*10\\^/'`
                        set tblBb     = \
                         `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f2 | sed -e 's/[eE]+*/\\*10\\^/'`
                        ## SATDW = BB*SATDK*(SATPSI/MAXSMC)
                        set newValue = \
                            `echo "scale=16; $tblBb*$tblSatDk*($tblSatPsi/$tblMaxSmc)" | bc`
                    endif 
                else 
                    set origValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f9`
                    set origValue = `printf '%.16f' $origValue`
                    set newValue = `echo "$origValue * $theSatDw" | bc`
                endif
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p SOILPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$9="  "newValue;print}' \
                    >> SOILPARM.TBL.NEW 
            end

            mv SOILPARM.TBL.NEW   SOILPARM.TBL
        endif ## 1D

        if ( $satDwDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 

    endif  # satDw




    #--------------------------------------------    
    #--------------------------------------------    
    ## VEGPARM.TBL

    ## rsMult ------------
    if ( `echo $assimOnlyVars | grep rsMult | wc -l` ) then 

        set rsMultDims = `ncks -m restart.assimOnly.nc | grep 'rsMult dimension' | wc -l`
        if ( $rsMultDims == 1 ) then
            set theRsMult = `ncdump -v rsMult restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`

            ## VEGPARM.TBL
            ## This is read in using commas. However, there were many commas missing. 
            ## Respecting the field width is for human readability. 
            set nlines = `cat VEGPARM.TBL | wc -l`
            \rm -rf VEGPARM.TBL.NEW
            cp VEGPARM.TBL VEGPARM.TBL.NEW

            ## these are checks that should be done to ensure proper formatting of 
            ## the vegparm table. THere were many missing commas, which will result in 
            ## issues, so I'm imposing a total comma count here.
            #for i in `seq 1 $nLinesVegparm`; do s#ed -n ${i}p VEGPARM.TBL | tr -C ',' ' ' | tr -d ' ' | wc -c; done
            set nCommas = `cat VEGPARM.TBL | tr -C ',' ' ' | tr -d ' ' | wc -c`
            if ( $nCommas != 2213 ) then
                echo '\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/'
                echo '!!!!!!!!!!!!    MALFORMED VEGPARM.TBL   !!!!!!!!!!!!   '
                echo '  OR it has changed since apply_assim.csh was written.'
                echo '!!!!!!!!!!!!    MALFORMED VEGPARM.TBL   !!!!!!!!!!!!   '
                echo '\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/'
                exit 86
            endif 

            # identify the lines with " RS " in them
            set whHeaders = `grep -n ' RS ' VEGPARM.TBL.NEW | cut -d':' -f1`

            foreach hh ( $whHeaders )
                set nLinesTable=`sed -n ${hh}p VEGPARM.TBL.NEW  | cut -d',' -f1`
                set whColRs=`sed -n ${hh}p VEGPARM.TBL.NEW | cut -d"'" -f2`
                set whColRs=`echo $whColRs | tr ' ' '\n' | grep -n RS | cut -d':' -f1`
                # since the dummy row numbers take up a column
                @ whColRs = $whColRs + 1

                foreach ll0 ( `seq 1 $nLinesTable` )
                    @ ll = $ll0 + $hh
                    set origValue = `sed -n ${ll}p VEGPARM.TBL.NEW | cut -d',' -f${whColRs}`
                    set origValue = `printf '%.16f' $origValue`
                    set newValue = `echo "$origValue * $theRsMult" | bc`
                    set newValue = `printf '%8f' $newValue`
                    # dealing with string literal is annoying. another reason to switch to bash?
                    set origLine = `sed -n ${ll}p VEGPARM.TBL.NEW | sed 's/ /\\ /g' | sed 's#/#\\/#g'`
                    set newLine = `echo $origLine | awk -v newValue=$newValue -v whColRs=$whColRs \
                      'BEGIN{FS=",";OFS=","};{$whColRs=newValue;print}'`
                    sed -i "${ll}s/.*/$newLine /" VEGPARM.TBL.NEW
                end 
            end

            mv VEGPARM.TBL.NEW VEGPARM.TBL

        endif ## 1D

        if ( $rsMultDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 

    #--------------------------------------------    
    #--------------------------------------------    
    ## MPTABLE.TBL

    ## ch2opMult ------------
    if ( `echo $assimOnlyVars | grep ch2opMult | wc -l` ) then 

        set ch2opMultDims = `ncks -m restart.assimOnly.nc | grep 'ch2opMult dimension' | wc -l`
        if ( $ch2opMultDims == 1 ) then
            set theCh2opMult = `ncdump -v ch2opMult restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`

            ## MPTABLE.TBL
            ## This is read in as a namelist, so repsecting the field widths 
            ## is only for human readability
            \rm -rf MPTABLE.TBL.NEW
            cp MPTABLE.TBL MPTABLE.TBL.NEW

            # identify the lines with " RS " in them
            set whCh2op = `grep -in 'ch2op' MPTABLE.TBL.NEW | cut -d':' -f1`

            foreach hh ( $whCh2op )
                set theLine = `sed -n ${hh}p MPTABLE.TBL.NEW | sed 's/ /\\ /g'`
                set theLine2 = `echo $theLine | cut -d'=' -f2`

                set outLine = `echo $theLine | cut -d'=' -f1`
                set outLine = "$outLine ="

                foreach ll0 ( `echo $theLine2 | tr ' ' 'j' | tr ',' ' ' | tr -d '\\,'` )
                    set fieldWidth = `echo $ll0 | wc -c`
                    @ fieldWidth = $fieldWidth - 1
                    set ll = `echo $ll0 | tr 'j' ' '`
                    set newValue = `echo "$ll * $theCh2opMult" | bc` 
                    set newValue = `printf "%${fieldWidth}f" $newValue | cut -c1-${fieldWidth}`
                    set outLine = `echo "${outLine}${newValue},"`
                    #echo ll, theCh2opMult, newValue: $ll, $theCh2opMult , $newValue
               end

                sed -i "${hh}s/.*/$outLine/" MPTABLE.TBL.NEW
            end

            mv MPTABLE.TBL.NEW MPTABLE.TBL

        endif ## 1D

        if ( $ch2opMultDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 


    #--------------------------------------------    
    #--------------------------------------------    
    ## GENPARM.TBL
    ## GENPARM must be below SOILPARM

    ## SLOPE -----------------
    if ( `echo $assimOnlyVars | grep slope | wc -l` ) then 
        set slopeDims = `ncks -m restart.assimOnly.nc | grep 'slope dimension' | wc -l`
        if ( $slopeDims == 1 ) then
            set theSlope = `ncdump -v slope restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            set theSlope = `echo $theSlope | cut -c1-5`
            set lineNumSlope = `grep -n SLOPE_DATA GENPARM.TBL | cut -d: -f1`
            @ setLineNum = $lineNumSlope + 2
            sed -i "${setLineNum}s/.*/$theSlope/" GENPARM.TBL
        endif 
        if ( $slopeDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif

    ## refDk ----------------- This has to be written AFTER SATDK
    ## if set to -1 in restart.assimOnly.nc, then use the satDk value in the SOILPARM.TBL
    if ( `echo $assimOnlyVars | grep refDk | wc -l` ) then 
        set refDkDims = `ncks -m restart.assimOnly.nc | grep 'refDk dimension' | wc -l`
        if ( $refDkDims == 1 ) then
            set theRefDk = `ncdump -v refDk restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## if refDk is -1 then we set refDk using the current SOILPARM.TBL value
            ## of satdk for silty clay loam
            if ( `echo "$theRefDk == -1" | bc` ) then
                set whSltClyLom=`grep -n 'SILTY CLAY LOAM' SOILPARM.TBL | cut -d: -f1`
                set whSltClyLom=$whSltClyLom[1]
                set theRefDk = `sed -n ${whSltClyLom}p SOILPARM.TBL | cut -d',' -f8`
            endif
            set lineNumRefDk = `grep -n REFDK_DATA GENPARM.TBL | cut -d: -f1`
            @ setLineNum = $lineNumRefDk + 1
            sed -i "${setLineNum}s/.*/$theRefDk/" GENPARM.TBL
        endif 
        if ( $refDkDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif

    ## refKdt ----------------- 
    if ( `echo $assimOnlyVars | grep refKdt | wc -l` ) then 
        set refKdtDims = `ncks -m restart.assimOnly.nc | grep 'refKdt dimension' | wc -l`
        if ( $refKdtDims == 1 ) then
            set theRefKdt = `ncdump -v refKdt restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            set lineNumRefKdt = `grep -n REFKDT_DATA GENPARM.TBL | cut -d: -f1`
            @ setLineNum = $lineNumRefKdt + 1
            sed -i "${setLineNum}s/.*/$theRefKdt/" GENPARM.TBL
        endif 
        if ( $refKdtDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif

    ## czil ------------
    if ( `echo $assimOnlyVars | grep czil | wc -l` ) then 

        set czilDims = `ncks -m restart.assimOnly.nc | grep 'czil dimension' | wc -l`
        if ( $czilDims != 1 ) then
            echo 'CZIL should be a scalar!'
            exit 3
        endif 

        set theCzil = `ncdump -v czil restart.assimOnly.nc | tail -n2 | head -1 | \
                         cut -d'=' -f2 | tr -d ' ;'`

        \rm -rf GENPARM.TBL.NEW
        cp GENPARM.TBL GENPARM.TBL.NEW
        # identify the lines with " RS " in them
        set whCzil = `grep -in 'czil_data' GENPARM.TBL.NEW | cut -d':' -f1`
        @ whCzilVal = $whCzil + 1
        sed -i "${whCzilVal}s/.*/$theCzil/" GENPARM.TBL.NEW

        mv GENPARM.TBL.NEW GENPARM.TBL
 
    endif  # czil

endif  # paramActive

exit 0

#!/bin/csh
@ errCode = 0
set theDate = `date +'%Y-%m-%d_%H:%M:%S_%Z'` || exit `expr $errCode +1`

## create the log file. 
set logFile = wrfHydroMetaLog.${theDate} || exit `expr $errCode +1`
rm -rf $logFile || exit `expr $errCode +1`
touch $logFile || exit `expr $errCode +1`

## log time 
## this first line/date in the file is included with each section header
## to help assist parsing the file.
echo $theDate >> $logFile || exit `expr $errCode +1`
echo >> $logFile

## log working dir
echo "run dir --- $theDate" >> $logFile
echo `pwd` >> $logFile || exit `expr $errCode +1`
echo  >> $logFile

## log namelist.hrldas
echo "namelist.hrldas --- $theDate" >> $logFile
cat namelist.hrldas >> $logFile || exit `expr $errCode +1`
echo >> $logFile

## log hydro.namelist
echo "hydro.namelist --- $theDate" >> $logFile
cat hydro.namelist >> $logFile || exit `expr $errCode +1`
echo >> $logFile

## log generic restart links
foreach rst (`ls restart*.nc`)
    if ( -e $rst ) then
	echo "$rst --- $theDate" >> $logFile
	readlink -e $rst >> $logFile
	echo >> $logFile
    endif
end

## log forcing link
set theInDir = `grep INDIR namelist.hrldas | egrep -v '^!' | tr -d " '" | tr -d '"' | cut -d'=' -f2 | cut -d'!' -f1`
echo "INDIR --- $theDate" >> $logFile
readlink -e $theInDir >> $logFile
echo >> $logFile

## log parameters
foreach TBL (`ls *.TBL`)
    echo "$TBL --- $theDate" >> $logFile
    cat $TBL >> $logFile
    echo >> $logFile
end

exit 0

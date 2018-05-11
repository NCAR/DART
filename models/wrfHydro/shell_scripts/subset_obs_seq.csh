#!/bin/csh

## This script make it easy to subset an obs_seq.out file based on a specified number of hours
## and the restart times in the current initialEnsemble/ folder (for the LSM component ens 1)
## It edits input.nml to do this. 
## The given input_obs_seq.out filename is appended with the date range in the file. 
## Todo: if the input/default input_obs_seq.out is a symlink: die or follow readlink?

set restartTime = `ncdump -v Times initialEnsemble/restart.0001.nc | tail -2 | head -1 | cut -d'"' -f2`
set restartTime2 = `echo $restartTime | cut -d: -f1 | tr '_' ' '`" hour"

set startTime = `date -ud "$restartTime2" +'%Y %m %d %H'`

if ($#argv != 1 & $#argv != 2) then
    echo "Usage: nHours [input_obs_seq.out]"
    echo "Start time is that of restart.0001.nc file in initialEnsemble:" $startTime
    echo "Arg 2 optional, defaults to the obs_seq.out file already set in the OBS_SEQUENCE_TOOL_NML in input.nml"
    exit 1
endif
echo "Start time is that of restart.0001.nc file in initialEnsemble:" $startTime

set nHours    = $1
if ($#argv == 2) set theObsSeq = $2

set endTime = `date -ud "$restartTime2 + $nHours hours" +"%Y %m %d %H"`

set startTimeCompact = `echo $startTime | tr -d ' '`
set endTimeCompact = `echo $endTime | tr -d ' '`

set gregStart = `./gregorian_time $startTime 00 00`
set gregStartDay = `echo $gregStart | cut -d' ' -f1`
set gregStartSec = `echo $gregStart | cut -d' ' -f2`

set gregEnd = `./gregorian_time $endTime 00 00`
set gregEndDay = `echo $gregEnd | cut -d' ' -f1`
set gregEndSec = `echo $gregEnd | cut -d' ' -f2`

@ obsSeqStartLine = `grep -n 'obs_sequence_tool_nml' input.nml | cut -d: -f1`

@ nInObsSeqLine = $obsSeqStartLine - 1 + \
      `tail -n+$obsSeqStartLine input.nml | head -13 | grep -n num_input_files | cut -d: -f1` 
@ inObsSeqLine = $obsSeqStartLine - 1 + \
      `tail -n+$obsSeqStartLine input.nml | head -13 | grep -n filename_seq | cut -d: -f1` 
@ outObsSeqLine = $obsSeqStartLine - 1 + \
      `tail -n+$obsSeqStartLine input.nml | head -13 | grep -n filename_out | cut -d: -f1` 

@ firstDayLine = $obsSeqStartLine - 1 + \
      `tail -n+$obsSeqStartLine input.nml | head -13 | grep -n first_obs_days | cut -d: -f1` 
@ firstSecLine = $obsSeqStartLine - 1 + \
      `tail -n+$obsSeqStartLine input.nml | head -13 | grep -n first_obs_seconds | cut -d: -f1` 

@ lastDayLine = $obsSeqStartLine - 1 + \
      `tail -n+$obsSeqStartLine input.nml | head -13 | grep -n last_obs_days | cut -d: -f1` 
@ lastSecLine = $obsSeqStartLine - 1 + \
      `tail -n+$obsSeqStartLine input.nml | head -13 | grep -n last_obs_seconds | cut -d: -f1` 

if (! $?theObsSeq) set theObsSeq = `tail -n+$inObsSeqLine input.nml | head -1 | cut -d'=' -f2 | tr -d "'" | tr -d '"'`
echo $theObsSeq

set nInObsSeq = `echo $theObsSeq | cut -d= -f2 | tr -d ' ' | tr ',' ' '`
set nInObsSeq = $#nInObsSeq

## sed -i unlinks the file replacing with a copy, so edit the target
set linkedFile = `readlink -e input.nml`
## the new file/ subset
set newObsSeqOutSub = "${theObsSeq}.${startTimeCompact}-${endTimeCompact}"
echo "Output file: ${theObsSeq}.${startTimeCompact}-${endTimeCompact}"

sed -i "${nInObsSeqLine}s/.*/   num_input_files    = 1/" $linkedFile
sed -i "${inObsSeqLine}s#.*#   filename_seq       = '${theObsSeq}'#" $linkedFile
sed -i "${outObsSeqLine}s#.*#   filename_out       = '${newObsSeqOutSub}'#" $linkedFile

sed -i "${firstDayLine}s/.*/   first_obs_days     = ${gregStartDay}/" $linkedFile
sed -i "${firstSecLine}s/.*/   first_obs_seconds  = ${gregStartSec}/" $linkedFile

sed -i "${lastDayLine}s/.*/   last_obs_days      = ${gregEndDay}/" $linkedFile
sed -i "${lastSecLine}s/.*/   last_obs_seconds   = ${gregEndSec}/" $linkedFile

./obs_sequence_tool

set failure = $?
if ( $failure ) exit $failure

ln -sf $newObsSeqOutSub obs_seq.out

echo "Output file: ${theObsSeq}.${startTimeCompact}-${endTimeCompact}"
echo "Is linked to obs_seq.out"

exit 0

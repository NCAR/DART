#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Since the files are already created with a 'time' record dimension,
# It is perhaps simpler to concatenate along the time dimension and
# then concatenate all the ensemble members. This unfortunately
# results in variables (member,time,<whatever>) rather than the desired
# (time,member,<whatever>), the dimensions must be permuted. Bummer.
# I did try to concatenate the members first (which creates a SECOND
# unlimited dimension) and it got really wierd trying to specify
# which variables get the member record dimension, which ones get the
# missing 'time' dimension, which dimension to concatenate along ...
# I am just living with the slight inefficiency. TJH

cat << EndOfFile >! time.template.cdl
netcdf time.template {
dimensions:
        time = UNLIMITED ; // (1 currently)
variables:
        int time(time) ;
		time:long_name = "valid output time" ;
		time:standard_name = "time" ;
		time:units = "minutes since 1970-01-01 00:00:00 UTC" ;
// global attributes:
                :version = "dummy for time" ;
data:
 time =  YYYYMMDD ;
}
EndOfFile

# requires 1 argument:
#    $1 - the ensemble member integer that is part of the directory name

if ($# != 1) then
   echo usage: $0 member_number
   exit -1
endif

set member = $1     # target date, YYYYMMDD

set MYEXT = `printf        %03d $member`
set MYDIR = `printf member_%03d $member`

cd ${MYDIR}

ln -s ../input.nml .

#==========================================================================
# 201306250100.CHRTOUT_DOMAIN1
#==========================================================================

set VARLIST = 'crs,Head,qBtmVertRunoff,qBucket,qSfcLatRunoff,q_lateral,streamflow,velocity'
\rm -f output.nc tempvar.nc

# add a time dimension to the variables that dont have one
foreach FILE ( *.CHRTOUT_DOMAIN1 )
   set TEMPFILE = ${FILE}.tempvar.nc
   ncecat -O -h -H -v ${VARLIST} -u time ${FILE} ${TEMPFILE}
   echo "copying $FILE to $TEMPFILE"
end

# concatenate input files without target_var
ls -1 *.CHRTOUT_DOMAIN1 | sort | ncrcat -h -H -x -v ${VARLIST} -O -o output.nc

# concatenate input files with only target_var
ls -1 *.CHRTOUT_DOMAIN1.tempvar.nc | sort | ncrcat -h -H -O -o tempvar.nc

# combine variables into output.nc, rename, and clean up.
ncks -A tempvar.nc output.nc
mv output.nc ../CHRTOUT_DOMAIN1.${MYEXT}.nc
\rm -f tempvar.nc *.tempvar.nc

#==========================================================================
# 201306250100.CHANOBS_DOMAIN1
#==========================================================================

set VARLIST = 'crs,streamflow'

# add a time dimension to the variables that dont have one
foreach FILE ( *.CHANOBS_DOMAIN1 )
   set TEMPFILE = ${FILE}.tempvar.nc
   ncecat -O -h -H -v ${VARLIST} -u time ${FILE} ${TEMPFILE}
   echo "copying $FILE to $TEMPFILE"
end

# concatenate input files without target_var
ls -1 *.CHANOBS_DOMAIN1 | sort | ncrcat -h -H -x -v ${VARLIST} -O -o output.nc

# concatenate input files with only target_var
ls -1 *.CHANOBS_DOMAIN1.tempvar.nc | sort | ncrcat -h -H -O -o tempvar.nc

# combine variables into output.nc, rename, and clean up.
ncks -A tempvar.nc output.nc
mv output.nc ../CHANOBS_DOMAIN1.${MYEXT}.nc
\rm -f tempvar.nc *.tempvar.nc

#==========================================================================
# 201306250100.LAKEOUT_DOMAIN1
#==========================================================================

set VARLIST = 'crs,water_sfc_elev,inflow,outflow'

# add a time dimension to the variables that dont have one
foreach FILE ( *.LAKEOUT_DOMAIN1 )
   set TEMPFILE = ${FILE}.tempvar.nc
   ncecat -O -h -H -v ${VARLIST} -u time ${FILE} ${TEMPFILE}
   echo "copying $FILE to $TEMPFILE"
end

# concatenate input files without target_var
ls -1 *.LAKEOUT_DOMAIN1 | sort | ncrcat -h -H -x -v ${VARLIST} -O -o output.nc

# concatenate input files with only target_var
ls -1 *.LAKEOUT_DOMAIN1.tempvar.nc | sort | ncrcat -h -H -O -o tempvar.nc

# combine variables into output.nc, rename, and clean up.
ncks -A tempvar.nc output.nc
mv output.nc ../LAKEOUT_DOMAIN1.${MYEXT}.nc
\rm -f tempvar.nc *.tempvar.nc

#==========================================================================
# HYDRO_RST.2013-06-01_01:00_DOMAIN1  no time dimension/vars in file
#==========================================================================

# add a time dimension and variable to each file
foreach FILE ( HYDRO_RST.*_DOMAIN1 )

   set TEMPFILE = ${FILE}.tempvar.nc
   ncecat -O -h -H -u time ${FILE} ${TEMPFILE}
   echo "copying ${MYEXT} ${FILE} to ${TEMPFILE}"

   # Create a time variable from the time in the filename
   set TIMESTRING = `echo $FILE:e | sed "s/_DOMAIN1//"`
   set OFFSET = `echo $TIMESTRING 0 -g | ../advance_time`
   set DAYS = $OFFSET[1]
   set SECS = $OFFSET[2]
   set YYYYMMDD = `echo "($DAYS*24*60) + ($SECS/60)" | bc`
   sed -e "s/YYYYMMDD/${YYYYMMDD}/g" ../time.template.cdl >! time.cdl
   ncgen -o mytime.nc time.cdl

   # finally append to temporary file
   ncks -h -H -A -v time mytime.nc ${TEMPFILE}
   \rm mytime.nc time.cdl
end

# this one is the wrong size
\rm -f HYDRO_RST.2013-06-01_00:00_DOMAIN1.tempvar.nc

# concatenate files
ls -1 HYDRO_RST.*.tempvar.nc | sort | ncrcat -h -H -O -o ../HYDRO_RST.${MYEXT}.nc

# clean up.
\rm -f *.tempvar.nc input.nml time.template.cdl

#==========================================================================

exit 0


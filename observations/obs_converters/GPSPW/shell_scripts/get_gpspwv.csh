#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Fetch hourly pwv data from the Suominet website and rename the files.
# Soyoung Ha (MMM/NCAR) Dec-12-2014

set  yr = 2012
set  mo = 05
set day = 121   # a reference day (e.g., total number of days up to 2012-05-01)
set idy = 146	# 2012-05-25
set edy = 152   # 2012-05-31

set region = conus	# conus or globe
set fhead = CsuPWVh
set f_new = GPSPW_conus
set subdr = ncConusHourly
set ftail = "00.0060"

set DIR_SUOMINET = http://www.suominet.ucar.edu/data
set DIR_LOCAL    = /glade/p/nmmm0024/syha/OBS_SEQ/GPSPW/data

setenv  WGET /usr/bin/wget

while ( $idy <= $edy )

  set dy = `expr $idy \- $day`
  set dd = `printf %02d $dy`

  foreach ih ( 00 06 12 18 )
  set fi = ${fhead}_${yr}.${idy}.${ih}.${ftail}_nc
  set fo = ${f_new}_${yr}${mo}${dd}${ih}.nc
  ${WGET} -nv -t 1 ${DIR_SUOMINET}/$subdr/$fi
  mv $fi ${DIR_LOCAL}/$fo
  ls -l ${DIR_LOCAL}/$fo
  end

@ idy++
end

exit 0



#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#----------------------------------------------------------------------

@ year = 1998
@ ens_member = 1
@ max_ensemble_mem = 80

set  INPUTDIR = /glade/scratch/thoar/POP11/6hourly
set OUTPUTDIR = /glade/proj3/DART/CAM_DATM/${year}

#----------------------------------------------------------------------

cd ${INPUTDIR}
mkdir -p ${OUTPUTDIR}

while ( ${ens_member} <= ${max_ensemble_mem} )

   set fbase = `printf CAM_DATM-%02d.cpl.ha2x1dx6h ${ens_member}`

   set   jan = `printf ${fbase}.%04d-01.nc ${year}`
   set   feb = `printf ${fbase}.%04d-02.nc ${year}`
   set   mar = `printf ${fbase}.%04d-03.nc ${year}`
   set   apr = `printf ${fbase}.%04d-04.nc ${year}`
   set   may = `printf ${fbase}.%04d-05.nc ${year}`
   set   jun = `printf ${fbase}.%04d-06.nc ${year}`
   set   jul = `printf ${fbase}.%04d-07.nc ${year}`
   set   aug = `printf ${fbase}.%04d-08.nc ${year}`
   set   sep = `printf ${fbase}.%04d-09.nc ${year}`
   set   oct = `printf ${fbase}.%04d-10.nc ${year}`
   set   nov = `printf ${fbase}.%04d-11.nc ${year}`
   set   dec = `printf ${fbase}.%04d-12.nc ${year}`

   set ofile = `printf ${fbase}.%04d.nc    ${year}`

   echo -n "Appending 6hourly ensemble member ${ens_member} ... "

   ncrcat -O ${jan} ${feb} ${mar} ${apr} ${may} ${jun} ${jul} ${aug} ${sep} ${oct} ${nov} ${dec} \
             ${OUTPUTDIR}/${ofile}  || exit ${ens_member}

   echo "done at "`date`

#  echo -n "Fixing the time:calendar attribute for member $ens_member ... "
#  ncatted -O -a calendar,time,o,c,gregorian $ofile  || exit $ens_member
#  echo "done at "`date`

   @ ens_member ++

end

exit 0



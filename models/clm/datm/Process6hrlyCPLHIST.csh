#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

set DATADIR = /glade/scratch/thoar/POP15/6hourly
set DATADIR = /glade/scratch/thoar/POP11/6hourly
cd ${DATADIR}

#----------------------------------------------------------------------
# Create a teeny netCDF file with the right time bounds.
# We start/stop CAM every 6 hours, and the coupler history
# files record the start time of the run as the
# 'units' portion of the "time" variable and the
# length of time 'into' a run as the value of "time". 
#
# Consequently, every time variable we generate is "0.25" 
# We would like to reference them to a common start.
#
# We will then overwrite the variables in the original
# netCDF files with the values in these variables.
#----------------------------------------------------------------------

@ year = 1999

@ start_month = 1
@ last_month = 12

# The proper offset for days in december is 334 for a noleap year.

set offset = ( 0 31 59 90 120 151 181 212 243 273 304 334 )

cat << ENDOFFILE >! bob.ncdump
netcdf bob {
dimensions:
        time = UNLIMITED ; // (1 currently)
        ntb = 2 ;
variables:
        double time(time) ;
                time:units = "days since ${year}-01-01 00:00:00" ;
                time:calendar = "noleap" ;
                time:bounds = "time_bnds" ;

        double time_bnds(time, ntb) ;

// global attributes:
                :file_version = "tim" ;
data:
 time = xxx.25, xxx.50, xxx.75, yyy.00;

 time_bnds =
  xxx.00, xxx.25,
  xxx.25, xxx.50,
  xxx.50, xxx.75,
  xxx.75, yyy.00;
}
ENDOFFILE

cat << LASTDAYOFYEAR >! lastday.ncdump
netcdf lastday {
dimensions:
        time = UNLIMITED ; // (1 currently)
        ntb = 2 ;
variables:
        double time(time) ;
                time:units = "days since ${year}-01-01 00:00:00" ;
                time:calendar = "noleap" ;
                time:bounds = "time_bnds" ;

        double time_bnds(time, ntb) ;

// global attributes:
                :file_version = "tim" ;
data:
 time = xxx.25, xxx.50, xxx.75, xxx.999;

 time_bnds =
  xxx.00, xxx.25,
  xxx.25, xxx.50,
  xxx.50, xxx.75,
  xxx.75, xxx.999;
}
LASTDAYOFYEAR

#----------------------------------------------------------------------
#
#----------------------------------------------------------------------

@ ens_member = 1

while ( $ens_member <= 80 )


   # single-POP15 for 2000, 2001, 2002
   set fbase = "FV2deg_Cplr_out_single-POP15-${ens_member}.cpl.ha2x1dx6h"

   # CAM_halo-O2-POP15 for 2003
   set fbase = "CAM_halo-O2-POP15-${ens_member}.cpl.ha2x1dx6h"

   # POP11 for 1997, 1998, 1999
   set fbase = "FV_greg_single-O2-POP11-${ens_member}.cpl.ha2x1dx6h"

   @ month = $start_month

   while ( $month <= $last_month )

      set ofile = `printf CAM_DATM-%02d.cpl.ha2x1dx6h.%04d-%02d.nc ${ens_member} ${year} ${month}`
      rm -f ${ofile}

      set myfiles = `printf ${fbase}.%04d-%02d ${year} ${month}`

      foreach FILE ( ${myfiles}-*.nc )

         #----------------------------------------------------------------
         # Extract the day from the YYYY-MM-DD.nc filename
         #----------------------------------------------------------------

         set FILEBASE = $FILE:r
         set MODEL_DATE_STR = `echo $FILEBASE:e | sed -e "s#-# #g"`
         set MODEL_DATE = `echo $MODEL_DATE_STR`
         @ day = `echo $MODEL_DATE[3] | sed -e "s#08#8#" -e "s#09#9#"`
         
         #----------------------------------------------------------------
         # Create a netCDF file with the right time attributes
         # The last day of the year must be 364.999 .NOT. 365.00
         # There are no leap years.
         #----------------------------------------------------------------

         @ time0 = $day + $offset[$month] - 1
         @ time1 = $time0 + 1

         if ( $time0 == 364 ) then
            sed -e "s#xxx#${time0}#g" lastday.ncdump >! new.foo
         else
            sed -e "s#xxx#${time0}#g" -e "s#yyy#${time1}#g" bob.ncdump >! new.foo
         endif

         ncgen -o bob.nc new.foo

         #----------------------------------------------------------------
         # unpack the netCDF file and replace the time variables
         #----------------------------------------------------------------

         set sequence = `printf ${fbase}.%04d-%02d-%02d.nc $year $month $day`

         echo "Replacing times in ${sequence}"

         if ( -e ${sequence}.gz ) gunzip ${sequence}.gz

         ncks -A -v time,time_bnds bob.nc ${sequence} || exit $day

      end

      set sequence = `printf ${fbase}.%04d-%02d $year $month`

      # Move any possible february 29th out of the way.
      if ( -e ${sequence}-29.nc && $month == 2 ) then
         mv ${sequence}-29.nc ${sequence}-29.nc.unwanted 
      endif

      # Combine all the timesteps into one
      echo -n "Combining into monthly file named ${ofile}"
      ncrcat -O ${sequence}-??.nc ${ofile} || exit ${ens_member}
      echo " done at "`date`

      # remove original datafile if space is at a premium
      # rm ${sequence}-??.nc

      @ month ++

   end

   @ ens_member ++

   rm -f temp.nc
end

rm -f bob.ncdump lastday.ncdump bob.nc foo new.foo

exit 0



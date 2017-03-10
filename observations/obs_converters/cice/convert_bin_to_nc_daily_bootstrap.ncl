; DART software - Copyright UCAR. This open source software is provided
; by UCAR, "as is", without charge, subject to all terms of use at
; http://www.image.ucar.edu/DAReS/DART/DART_download
;
; DART $Id$

; this script reads in the daily fields of the bootstrap data one year at a time and 
; archive as netcdf files with metadata

;============data description from CC matlab script=================================
; Data are stored as flat two-byte integers representing sea ice concentration
; values. The sea ice concentration data values are packed into integer
; format by multiplying the original sea ice concentration values by
; 10. These values range from 0 to 1,000, with land registered as -800 and
; -100 being the northern hemisphere hole poleward of 85 degrees for the
; SMMR data and poleward of 87 degrees for the SSM/I data, which the
; satellite cannot cover due to orbit inclination.
; little-endian byte order no mention of hearder
; One-byte (scaled, unsigned integer) flat binary arrays preceded by a 300-byte header;
; *uint8
;=====================================================================================

;obs located in /home/disk/oly26/bitz/observations/seaice/$year
;

nlat = 448
nlon = 304
; longitude ranges from -180 to 180

setfileoption("bin","ReadByteOrder","LIttleEndian")
lat = fbindirread("/home/disk/eos11/bitz/observations/seaice/data/psn25lats_v3.dat",0,(/nlat,nlon/),"integer")/100000.
lon = fbindirread("/home/disk/eos11/bitz/observations/seaice/data/psn25lons_v3.dat",0,(/nlat,nlon/),"integer")/100000.

lat!0 = "ny"
lat!1 = "nx"

lon!0 = "ny"
lon!1 = "nx"

lat@units = "degrees_north"
;lon@units = "degrees_east"
printVarSummary(lon)
printVarSummary(lat)

years = 1990
yeare = 1991
nyear = yeare-years + 1

dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
doml= (/31,29,31,30,31,30,31,31,30,31,30,31/)

filedir = "/home/disk/oly26/bitz/observations/seaice/"
outdir  = "/home/disk/eos8/yz4362/Work/observations/seaice/daily/"

do iyear = 0, nyear -1
  
   year = years + iyear

   do imonth = 0, 11

       month = sprinti("%0.2i",imonth+1)

       days=0
    
       if (year%4 .eq.0)then
  
           data = new((/doml(imonth),nlat,nlon/),float)
           daye = doml(imonth) - 1   
           time = ispan(0,daye,1)
  
       else
  
           data = new((/dom(imonth),nlat,nlon/),float)
           daye = dom(imonth) -1
           time = ispan(0,daye,1)

       end if

       time!0 = "time"
       time&time = time
       time@unit = "days since "+year+"-"+month+"-01"

       do iday =days,daye

          day = sprinti("%02i",iday + 1)

          filename = systemfunc("ls "+filedir +year+"/bt_"+year+month+day+"_*_v02_n.bin")

          setfileoption("bin","ReadByteOrder","LittleEndian")         
          temp = fbindirread(filename,0,-1,"short")
          temp_data= onedtond(temp,(/nlat,nlon/))

          data(iday,:,:) = temp_data/10.0

          delete(temp)
          delete(temp_data)
      
       end do ; end of the day loop

       data = where(data.gt.100.0,data@_FillValue,data)

       data = where(data.lt.0.0, data@_FillValue,data)

       printVarSummary(time)
       printVarSummary(data)
       data!0 = "time"
       data&time = time
       data!1 = "ny"
       data!2 = "nx"

   
       data@units  = "percentage"
       data@longname = "sea ice concentration"
 
       fileout = outdir + "bt_"+year+month+"_v02_n.nc"
       system("rm -rf "+fileout)
       ncdfout = addfile(fileout,"c")

       setfileoption(ncdfout,"DefineMode",True)

       fAtt = True
       fAtt@title        ="SSMR SSM/I sea ice concentration using bt algorithm"
       fAtt@source_file  =" sea ice concentration binary data; daily"
       fAtt@Conventions  ="None"
       fAtt@creation_data=systemfunc("date")

       dimNames    =(/"time","ny","nx"/)
       dimSizes    =(/-1,nlat,nlon/)
       dimUnlim    =(/True,   False,  False/)
       filedimdef(ncdfout,dimNames,dimSizes,dimUnlim)

       filevardef(ncdfout, "time", typeof(time), getvardims(time))
       filevardef(ncdfout, "lat",  typeof(lat),  getvardims(lat))
       filevardef(ncdfout, "lon",  typeof(lon),  getvardims(lon))
       filevardef(ncdfout, "sic",  typeof(data), getvardims(data))

       filevarattdef(ncdfout,"time",time)
       filevarattdef(ncdfout,"lat" ,lat)
       filevarattdef(ncdfout,"lon" ,lon)
       filevarattdef(ncdfout,"sic" ,data)

       setfileoption(ncdfout,"DefineMode",False)

       ncdfout->time = (/time/)
       ncdfout->lat  = (/lat/)
       ncdfout->lon  = (/lon/)
       ncdfout->sic  = (/data/)

       printVarSummary(data)
       delete(data)
       delete(time)

     end do   ; end of the month loop

 end do     ; end of the year loop

; <next few lines under version control, do not edit>
; $URL$
; $Revision$
; $Date$

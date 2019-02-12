; DART software - Copyright UCAR. This open source software is provided
; by UCAR, "as is", without charge, subject to all terms of use at
; http://www.image.ucar.edu/DAReS/DART/DART_download
;
; DART $Id$

case = "free_rsnw_2005to2010"

tmember = 30 ;pick one of the ensemble members as truth

dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)

year = 2005

filedir = "/glade/scratch_new/yfzhang/archive/cesm2_0_alpha06n/truth/"

outdir  = "/$WORK/observations/syn/cesm2/"+case+"/t"+tmember+"/"

ncdfhead = addfile("/$WORK/headfiles/cice_hist.nc","r")

lat = ncdfhead->TLAT
lon = ncdfhead->TLON

printVarSummary(lat)

ni = 320 
nj = 384

lat!0 = "nj"
lat!1 = "ni"

lon!0 = "nj"
lon!1 = "ni"

ncat = 5

miu = 0
sigma = 0.1

rseed1 = toint(systemfunc(" date +%s"))
rseed2 = toint((12345l*rseed1)%2147483398l)
random_setallseed(rseed1, rseed2)


do imonth =  3,9

   month = sprinti("%0.2i",imonth+1)
   
   if (year.eq.2005 .and. imonth.eq.3) then
      days = 1
   else
      days = 0
   end if

   do iday = days, dom(imonth)-1

      day  = sprinti("%0.2i",iday+1)

      date = year+"-"+month+"-"+day

      ncdf_t = addfile(filedir+case+".cice_00"+tmember+".r."+date+"-00000.nc","r")

      aicen  = doubletofloat(ncdf_t->aicen)

      aice   = dim_sum_n(aicen,0)

      aice_obs = new((/nj,ni/),float)

      do ilat = 0, nj -1
            
          do ilon = 0, ni - 1

                 noise_add = new(1,float)

		 if (aice(ilat,ilon).gt.0.01) then

                 if (.not.ismissing(lat(ilat,ilon)) .and. lat(ilat,ilon).gt.40) then

                    if (aice(ilat,ilon).gt.0) then
 
		     err = aice(ilat,ilon)*sigma
     
		     noise_add = random_normal(miu,err,1)

                     aice_obs(ilat,ilon) = noise_add + aice(ilat,ilon)

                    else

                     aice_obs(ilat,ilon) = aice(ilat,ilon)

                    end if

                 end if
                 end if

           end do

       end do


       aice_obs!0 = "nj"
       aice_obs!1 = "ni"
        
       aice_obs = where(aice_obs.gt.1,1,aice_obs)
       aice_obs = where(aice_obs.lt.0,0,aice_obs)


       fileout = outdir+"/synthetic_obs.aice."+date+".nc"

        system("rm -rf "+fileout)

        ncdfout = addfile(fileout,"c")

	fAtt = True
	fAtt@title        ="Obs for perfect model experiments; only observations northern to 40 degreen North are converted"
	fAtt@source_file  ="Restart files from a CICE5 openloop case"
	fAtt@Conventions  ="None"
	fAtt@creation_data=systemfunc("date")

	dimNames    =(/"nj","ni"/)
	dimSizes    =(/nj,ni/)
	dimUnlim    =(/ False,  False/)
	filedimdef(ncdfout,dimNames,dimSizes,dimUnlim)

	filevardef(ncdfout, "lat",  typeof(lat),  getvardims(lat))
	filevardef(ncdfout, "lon",  typeof(lon),  getvardims(lon))
	filevardef(ncdfout, "aice",  typeof(aice_obs), getvardims(aice_obs))

	filevarattdef(ncdfout,"lat" ,lat)
	filevarattdef(ncdfout,"lon" ,lon)
	filevarattdef(ncdfout,"aice" ,aice_obs)

	setfileoption(ncdfout,"DefineMode",False)

	ncdfout->lat  = (/lat/)
	ncdfout->lon  = (/lon/)
	ncdfout->aice  = (/aice_obs/)
	delete(aice_obs)

   end do

   print(year + month+" finished")

end do

; <next few lines under version control, do not edit>
; $URL$
; $Revision$
; $Date$


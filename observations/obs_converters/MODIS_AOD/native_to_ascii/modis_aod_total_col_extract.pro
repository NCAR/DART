;=======================================================================
; subroutines for modis_extract_hdf
;=======================================================================
; Get the Vdata information from an HDF file
function get_vd, filename, varname
   file_id = hdf_open(filename, /read)
   vd_ref = hdf_vd_find(file_id, strtrim(varname,2))
   vdata=hdf_vd_attach(file_id,vd_ref)
   nread=hdf_vd_read(vdata,var)
   hdf_vd_detach, vdata
   return, var
end
;
; Get the Scientific Data information from an HDF file
function read_modis, filename, varname
   sd_id = hdf_sd_start(filename,/read)
   index = hdf_sd_nametoindex(sd_id, varname)
   sds_id = hdf_sd_select(sd_id,index)
   hdf_sd_getdata, sds_id, dat
   hdf_sd_endaccess, sds_id
   hdf_sd_end, sd_id
   return, dat
end
;
; Get the date and time
function date_time,year,month,day,hour,inc,yy_out,mm_out,dd_out,hh_out
  if( (year mod 4) eq 0 ) then begin
     month_array=[31,29,31,30,31,30,31,31,30,31,30,31]
  endif else begin
     month_array=[31,28,31,30,31,30,31,31,30,31,30,31]
  endelse
  yy_out=year
  mm_out=month
  dd_out=day
  hh_out=hour+inc
  if (hh_out lt 0) then begin
     hh_out=24+hh_out
     dd_out=day-1
  endif
  if (dd_out eq 0) then begin
     mm_out=mm_out-1
     if (mm_out ne 0) then begin
        dd_out=month_array(mm_out-1)
     endif else begin
        yy_out=year-1
        mm_out=12
        dd_out=31
     endelse
  endif
  rc=0
  return, rc
end         
;     
;=======================================================================
; main IDL routine
; needs read_modis, get_vd, 
;=======================================================================
pro modis_aod_total_col_extract, dir_in, outf, year, month, day, $
hour, window, lon_min, lon_max, lat_min, lat_max 
;=======================================================================
!EXCEPT = 0
;
   print, 'IDL dir in   ', dir_in
   print, 'IDL file out ', outf
   print, 'IDL year     ', year
   print, 'IDL month    ', month
   print, 'IDL day      ', day
   print, 'IDL hour     ', hour
   print, 'IDL window   ', window
   print, 'IDL lon min  ', lon_min
   print, 'IDL lon max  ', lon_max
   print, 'IDL lat min  ', lat_min
   print, 'IDL lat max  ', lat_max
;
; QA/QC settings
   aod_min=0.055
   aod_max=1.500
   cld_qa_crit=0
   lnd_qa_crit=1
   ocn_qa_crit=1
   cld_qa=0
   lnd_qa=1
   ocn_qa=1
   debug = 0
;
; Define MODIS input file name
   modis_input_dir=dir_in
   modis_output_file=outf
   openw, unit, modis_output_file, /get_lun
;
; Print file information
   print, 'IDL - MODIS input directry: ', modis_input_dir
   print, 'IDL - MODIS output ascii file: ', modis_output_file
;
; Set days-of-the-month array
   for iwin=-window, window do begin
      rc=date_time(year,month,day,hour,iwin,yy_out,mm_out,dd_out,hh_out)
      print, iwin,yy_out,mm_out,dd_out,hh_out
      if( (yy_out mod 4) eq 0. ) then begin
         month_array=[31,29,31,30,31,30,31,31,30,31,30,31]
         cumdays=[31,60,91,121,152,182,213,244,274,305,335,366]
      endif else begin
         month_array=[31,28,31,30,31,30,31,31,30,31,30,31]
         cumdays=[31,59,90,120,151,181,212,243,273,304,334,365]
      endelse
;
; Calculate time factor
      time_factor = 0
      for kyear=1993, yy_out-1 do begin
         ireg=1
         ileap=0
         if((kyear mod 4) eq 0) then begin 
            ireg=0
            ileap=1
         endif
         time_factor=time_factor+ireg*365+ileap*366
      endfor
;
      act_cum=dd_out
      if (mm_out-2 ge 0) then begin
          act_cum=cumdays(mm_out-2)+dd_out
      endif
      mm_out_str=string(mm_out)
      if(mm_out lt 10) then begin
         mm_out_str='0'+strcompress(string(mm_out),/REMOVE_ALL)
      endif
      dd_out_str=string(dd_out)
      if(dd_out lt 10) then begin
         dd_out_str='0'+strcompress(string(dd_out),/REMOVE_ALL)
      endif
      hh_out_str=string(hh_out)
      if(hh_out eq 0) then begin
         hh_out_str='00'
      endif else begin
         if(hh_out lt 10) then begin
            hh_out_str='0'+strcompress(string(hh_out),/REMOVE_ALL)
         endif      
      endelse
      act_cum_str=string(act_cum)
      if(act_cum lt  100) then begin
         act_cum_str='0'+strcompress(string(act_cum),/REMOVE_ALL)
      endif
      if(act_cum lt 10) then begin
         act_cum_str='00'+strcompress(string(act_cum),/REMOVE_ALL)
      endif
;     
      path_dir1=modis_input_dir + '/' + strcompress(string(yy_out)) + mm_out_str + $
      dd_out_str + '00/'
      path_file1='MYD04_L2.A' + strcompress(string(yy_out)) + act_cum_str + $
      '.' + hh_out_str + '*'
      path_file2='MOD04_L2.A' + strcompress(string(yy_out)) + act_cum_str + $
      '.' + hh_out_str + '*'
      path_dir=strcompress(path_dir1, /REMOVE_ALL)
      path_file11=strcompress(path_file1, /REMOVE_ALL)
      path_file22=strcompress(path_file2, /REMOVE_ALL)
;
      filelist1=file_search(path_dir,path_file11)
      nfiles1=size(filelist1)
      filelist2=file_search(path_dir,path_file22)
      nfiles2=size(filelist2)
      if(nfiles1(0) eq 0) then begin
         filelist=filelist2
      endif
      if(nfiles2(0) eq 0) then begin
         filelist=filelist1
      endif
      if(nfiles1(0) ne 0 && nfiles2(0) ne 0) then begin
         filelist=[filelist1, filelist2]
      endif
      nfiles=size(filelist)
      nlimit=-1
      if (nfiles(0) gt 0.) then begin
         nlimit=nfiles(1)-1
      endif
      scale_factor=0.001
      sec2day=1./(24.*60.*60.)
      otype=1
      vert=0
;
      for ifile=0, nlimit do begin
         modis_input_file=strcompress(filelist(ifile),/REMOVE_ALL)
         print, 'file ',strcompress(modis_input_file)
;
; Read Scan Start Time
         name='Scan_Start_Time'
         time=read_modis(modis_input_file, name)
;
; Read Latitude;
         name = 'Latitude'
         lat=read_modis(modis_input_file, name)
;
; Read Longitude;
         name='Longitude'
         lon=read_modis(modis_input_file, name)
;
; Read Cloud Mask QA;
;         name='Cloud_Mask_QA'
;         cldmask=read_modis(modis_input_file, name)
;          modis_atmos, modis_input_file, name, cldmask
;
; Read Quality Assurance Land;
         name='Quality_Assurance_Land'
         lndmask=read_modis(modis_input_file, name)
;
; Read Quality Assurance Ocean;
         name='Quality_Assurance_Ocean'
         ocnmask=read_modis(modis_input_file, name)
;
; Read Optical Depth Land And Ocean;
         name='Optical_Depth_Land_And_Ocean'
         tau=read_modis(modis_input_file, name)*scale_factor
;
; Read Deep Blue Aerosol Optical Depth 550 Land;
         name='Deep_Blue_Aerosol_Optical_Depth_550_Land'
         taud=read_modis(modis_input_file, name)*scale_factor
;
; Get size of the MODIS input arrays
         dims=size(time)
         for i=0, dims(1)-1 do begin
            for j=0, dims(2)-1 do begin
;
; unpack cld_qa
;               bin=dec2bin(cldmask(i,j),8)
;               mask=bin2dec(num2str(bin(8)))
;               cld_qa=bin2dec(num2str(bin([6:7])))
;               surf_flag=bin2dec( num2str(bin([2:3])))           
;               nbyte=0
;               nbit1=1
;               nbit2=2
;               cld_qa=ishft(cldmask(i,j,nbyte) and 2^nbit1,-nbit1)
;               if(nbit2 gt 0.) then cld_qa=ishft(cldmask(i,j,nbyte) and 2^nbit1+2^nbit2,-nbit1)
;
; unpack qa_lnd
;               bin=dec2bin(lndmask(i,j,5),8)
;               dbrc=bin2dec(num2str(bin([2:3])));               
;               dbat=bin2dec(num2str(bin([4:5])));
;               dbcf=bin2dec(num2str(bin([6:7])));
;               dbuf=bin2dec(num2str(bin([8])));                               
;
; unpack qa_ocn
;               bin=dec2bin(ocnmask(i,j,5),8)
;               dbrc=bin2dec(num2str(bin([2:3])));               
;               dbat=bin2dec(num2str(bin([4:5])));
;               dbcf=bin2dec(num2str(bin([6:7])));
;               dbuf=bin2dec(num2str(bin([8])));                               
;
               aod=tau(i,j)
               aodd=taud(i,j)
               aodlon=lon(i,j)
               aodlat=lat(i,j)
               aodtim=time(i,j)
;
               tmp=(time(i,j)-float(time_factor)/sec2day)*sec2day + 1.
               for itim=0,11 do begin
                  if(float(cumdays(itim)) ge tmp) then break
               endfor      
               month_ot=itim+1
               if (month_ot eq 1) then begin
                  days=tmp
               endif else begin
                  days=(tmp-float(cumdays(month_ot-2)))
               endelse
               day_ot=floor(days)
               hours=(days-float(day_ot))*24. 
               hour_ot=floor(hours)
               minutes=(hours-float(hour_ot))*60.
               minute=floor(minutes)
               seconds=(minutes-float(minute))*60.
               second=floor(seconds)
;
; Apply obs window test
               if ((iwin eq -window && hour_ot ge hh_out) || $
               (iwin gt -window && iwin lt window) || $
               (iwin eq window && hour_ot lt hh_out)) then begin
;                     if(aod gt 0) then begin
;                        print, 'DATA PRINT',unit, otype, aodlat, aodlon, vert, $ 
;                        yy_out, month_ot, day_ot, hour_ot, minute, second, aod, aodd
;                     endif
;
; Apply qa and domain test
                  if ( aod gt aod_min && aod lt aod_max && aodtim gt 0 && $ 
                  aodlon ge lon_min && aodlon le lon_max && $
                  aodlat ge lat_min && aodlat le lat_max && $
                  cld_qa eq cld_qa_crit ) then begin
;
; AOD uncertainty assignment
                     if (aod lt 0.2) then begin
                        aoderr = 0.1;
                     endif
                     if (aod ge 0.2 && aod le 1.4) then begin
                           aoderr=0.05+0.20*aod
                     endif
                     if (aod gt 1.4) then begin
                        aoderr=0.2+0.4*aod
                     endif
;
; Write AOD data to output file
                     print, 'SAVED OBS ',otype, aodlat, aodlon, vert, $ 
                     yy_out, month_ot, day_ot, hour_ot, minute, second, aod, $
                     aoderr, format='(i2,2f10.4,7i5,2f8.4)'
                     printf, unit, otype, aodlat, aodlon, vert, $ 
                     yy_out, month_ot, day_ot, hour_ot, minute, second, aod, $
                     aoderr, format='(i2,2f10.4,7i5,2f8.4)'
                  endif
               endif
            endfor
         endfor
      endfor
   endfor
   close, unit
end              

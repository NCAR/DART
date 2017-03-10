; convert_hdf_to_netcdf.pro
;
; DART $Id$
;
; PURPOSE:
;   A simple reader for EOS-Aura Microwave Limb Sounder (MLS) Level 2 
;   Geophysical Product (L2GP version 2).
;   
;
;
; AUTHOR:
;   Modified by N Pedatella to output the data into a netcdf file
;     calling sequence is convert_hdf_to_netcdf, "filename",data,"outputfile"
;
;
;
;   Young-In Won [04/18/2008]
;   Goddard Earth Science Data Information Service Center (GES DISC)
;   NASA/Goddard Space Flight Center
;   Modified from ReaL2GP_STD.pro (from JPL)
;
; CALLING SEQUENCE:
;   IDL> read_mls_l2.pro, "filename", data
;    Example, IDL> readl2, "MLS-Aura_L2GP-BrO_v02-21-c01_2006d054.he5", data
;
;   IDL> help, data, /struct 		; shows the structure of data
;   IDL> print, data.ntimes 		; prints total profile number
;   IDL> print, data.time 		; prints time for each profile (TAI time)
;   IDL> print, data.nlevels 		; prints number of pressure levels
;   IDL> print, data.pressure 	; prints pressure level
;   IDL> print, data.longitude 	; prints longitude of each profile
;   IDL> print, data.l2gpvalue 	; prints geophysical values for each profile
;    ... so on
;
; REQUIRES:
;   IDL 6.1 or Greater
;
; INPUT PARAMETERS:
;   filename      - A scalar string that indicates the MLS L2GP file
;                   to be read (in HDF-EOS 5)
;   data	  - constructed data to be returned
;
; RETURN VALUE:
;   a data structure contains the following fields:
;   swathName, number of sample times (nTimes), number of vertical levels (nLevels), 
;   pressure, latitude, longitude, time, localSolarTime, solarZenithAngle, 
;   lineOfSightAngle, orbitGeodeticAngle, chunkNumber, l2gpValue, l2gpPrecision, 
;   status, quality, attributes.  If the any of these values are not found in the
;   file, they are replaced by a scalar value of the proper datatype.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 Pro convert_hdf_to_netcdf, filename, result, outfile

  COMPILE_OPT IDL2

  ;; Open hdf5 file and get fileID
  fileID = H5F_Open(filename)
  ;; Retreive the swath name(s). If more than one swath, select which swath to choose
  groupName = 'HDFEOS/SWATHS'
  print, groupName
  noSwaths = H5G_Get_NMembers(fileID, groupName)
  if (noswaths gt 1) then begin
  for swathindex=0, noswaths-1 do begin
  thisName = H5G_Get_Member_Name(fileID, groupName, swathIndex)
  print,strtrim(swathindex, 2), ') Swath= ', thisname
  endfor
  ;read, 'Which swath do you want? ', ans_swath
ans_swath = 'Temperature'  
swathname=H5G_Get_Member_Name(fileID, groupName, ans_swath)
  endif else begin
  swathname=H5G_Get_Member_Name(fileID, groupName, 0)
  endelse
  ;
  ; Get the geolocation fields info (Time, Latitude, Longitude, Pressure, etc.)
  ;
  path2=groupname+'/'+swathname+'/Geolocation Fields'
  geo_members=H5G_GET_NMEMBERS(fileid,path2)
  geo=strarr(geo_members)
  for jj= 0, geo_members-1 do begin
  geo[jj]=H5G_Get_Member_Name(fileid, path2, jj)
  endfor
  
  ; Reading geolocation Fields
  path3=groupname+'/'+swathname+'/Geolocation Fields/'+'Time'
  dsID = H5D_Open(fileid, path3)
  time = H5D_Read(dsiD)
  ntimes = N_Elements(time)
  H5D_Close, dsID
  path3=groupname+'/'+swathname+'/Geolocation Fields/'+'Latitude'
  dsID = H5D_Open(fileid, path3)
  latitude = H5D_Read(dsiD)
  H5D_Close, dsID
  path3=groupname+'/'+swathname+'/Geolocation Fields/'+'Longitude'
  dsID = H5D_Open(fileid, path3)
  longitude = H5D_Read(dsiD)
  H5D_Close, dsID
  
  ; For some swaths (e.g. O3 column-GEOS5), there is no levels
  dummy = Where(geo EQ 'Pressure', cnt)
  IF cnt GT 0 THEN BEGIN
    path3=groupname+'/'+swathname+'/Geolocation Fields/'+'Pressure'
    dsID = H5D_Open(fileid, path3)
    pressure = H5D_Read(dsID)
    H5D_Close, dsID
    nLevels = N_Elements(pressure)
  ENDIF ELSE BEGIN
    nLevels = 1
    pressure = [0.0]
  ENDELSE

  path3=groupname+'/'+swathname+'/Geolocation Fields/'+'LocalSolarTime'
  dsID = H5D_Open(fileid, path3)
  localsolartime = H5D_Read(dsiD)
  H5D_Close, dsID
  path3=groupname+'/'+swathname+'/Geolocation Fields/'+'SolarZenithAngle'
  dsID = H5D_Open(fileid, path3)
  solarzenithangle = H5D_Read(dsiD)
  H5D_Close, dsID
  path3=groupname+'/'+swathname+'/Geolocation Fields/'+'LineOfSightAngle'
  dsID = H5D_Open(fileid, path3)
  lineofsightangle = H5D_Read(dsiD)
  H5D_Close, dsID
  path3=groupname+'/'+swathname+'/Geolocation Fields/'+'OrbitGeodeticAngle'
  dsID = H5D_Open(fileid, path3)
  orbitgeodeticangle = H5D_Read(dsiD)
  H5D_Close, dsID
  path3=groupname+'/'+swathname+'/Geolocation Fields/'+'ChunkNumber'
  dsID = H5D_Open(fileid, path3)
  ChunkNumber = H5D_Read(dsiD)
  H5D_Close, dsID
 
  ;
  ; Get the data fields info (L2geophysicla Value, precision, status, quality, convergence)
  ;
  path4=groupname+'/'+swathname+'/Data Fields'
  data_members=H5G_GET_NMEMBERS(fileid,path4)
  for kk= 0, data_members-1 do begin
  dat=H5G_Get_Member_Name(fileid, path4, kk)
  endfor
  
  ;Reading data Fields (and attributes for the geophysical value only)
  atts = ''
  path3=groupname+'/'+swathname+'/Data Fields/'+'L2gpValue'
  dsID = H5D_Open(fileid, path3)
  l2gpvalue = H5D_Read(dsiD)
    FOR i = 0, H5A_Get_Num_Attrs(dsId) - 1 DO BEGIN
    attId = H5A_Open_Idx(dsId, i)
    atts = i EQ 0 ? Create_Struct(H5A_Get_Name(attId), (H5A_Read(attId))[0]) : $
                    Create_Struct(atts, H5A_Get_Name(attId), (H5A_Read(attId))[0])
    H5A_Close, attId
  ENDFOR
  H5D_Close, dsID
  
  path3=groupname+'/'+swathname+'/Data Fields/'+'L2gpPrecision'
  dsID = H5D_Open(fileid, path3)
  L2gpPrecision = H5D_Read(dsiD)
  H5D_Close, dsID
  
  path3=groupname+'/'+swathname+'/Data Fields/'+'Status'
  dsID = H5D_Open(fileid, path3)
  status = H5D_Read(dsiD)
  H5D_Close, dsID
  
  path3=groupname+'/'+swathname+'/Data Fields/'+'Quality'
  dsID = H5D_Open(fileid, path3)
  quality = H5D_Read(dsiD)  
  nLevels = N_Elements(pressure)
  H5D_Close, dsID
  
  path3=groupname+'/'+swathname+'/Data Fields/'+'Convergence'
  dsID = H5D_Open(fileid, path3)
  convergence = H5D_Read(dsiD)  
  H5D_Close, dsID

  ; Close file
  H5F_Close, fileID


  ; Construct the result
  result = {swathName:swathName, nTimes:nTimes, nLevels:nLevels, $
            pressure:pressure,latitude:latitude,longitude:longitude, $       
            time:time,localSolarTime:localSolarTime, $    
            solarZenithAngle:solarZenithAngle, $
            lineOfSightAngle:lineOfSightAngle, $      
            orbitGeodeticAngle:orbitGeodeticAngle, $    
            chunkNumber:chunkNumber, $
            l2gpValue:l2gpValue, $
            l2gpPrecision:l2gpPrecision, $
            status:status, $
            quality:quality, $
            convergence:convergence, $
            attributes:atts $
           }


;   NMP - output the data into a netcdf file
    id = NCDF_CREATE(outfile, /CLOBBER) ; open the netcdf file
    t_dimid = NCDF_DIMDEF(id,'times',ntimes)
    l_dimid = NCDF_DIMDEF(id,'levels',nlevels)
    t_vid = NCDF_VARDEF(id,'Time',[t_dimid],/DOUBLE)
    l_vid = NCDF_VARDEF(id,'Pressure',l_dimid,/FLOAT)
    lat_vid = NCDF_VARDEF(id,'Latitude',t_dimid,/FLOAT)
    lon_vid = NCDF_VARDEF(id,'Longitude',t_dimid,/FLOAT)
    lst_vid = NCDF_VARDEF(id,'LocalSolarTime',t_dimid,/FLOAT)
    temp_vid = NCDF_VARDEF(id,'Temperature',[l_dimid,t_dimid],/FLOAT)
    prec_vid = NCDF_VARDEF(id,'Precision',[l_dimid,t_dimid],/FLOAT)
    qual_vid = NCDF_VARDEF(id,'Quality',t_dimid,/FLOAT)
    con_vid = NCDF_VARDEF(id,'Convergence',t_dimid,/FLOAT)
    NCDF_CONTROL, id, /ENDEF ; TAKE OUT OF DEF MODE
    NCDF_VARPUT, id,t_vid,time
    NCDF_VARPUT, id,l_vid,pressure
    NCDF_VARPUT, id,lat_vid,latitude
    NCDF_VARPUT, id,lon_vid,longitude
    NCDF_VARPUT, id,lst_vid,localsolartime
    NCDF_VARPUT, id,temp_vid,l2gpvalue
    NCDF_VARPUT, id, prec_vid, l2gpPrecision
    NCDF_VARPUT, id, qual_vid, quality
    NCDF_VARPUT, id, con_vid, convergence
    NCDF_CLOSE, id ;close the file

END

; <next few lines under version control, do not edit>
; $URL$
; $Revision$
; $Date$

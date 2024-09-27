function omi_hcho_total_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
%
% Get file list and number of files
  wyr_mn=str2double(cwyr_mn);
   wmn_mn=str2double(cwmn_mn);
   wdy_mn=str2double(cwdy_mn);
   whh_mn=str2double(cwhh_mn);
   wmm_mn=str2double(cwmm_mn);
   wss_mn=str2double(cwss_mn);
   wyr_mx=str2double(cwyr_mx);
   wmn_mx=str2double(cwmn_mx);
   wdy_mx=str2double(cwdy_mx);
   whh_mx=str2double(cwhh_mx);
   wmm_mx=str2double(cwmm_mx);
   wss_mx=str2double(cwss_mx);
   nx_mdl=str2double(cnx_mdl);
   ny_mdl=str2double(cny_mdl);
%
   command=strcat('rm'," ",'-rf'," ",fileout);
   [status]=system(command);
   fid=fopen(fileout,'w');
%
   command=strcat('ls'," ",'-1'," ",filein,'*');
   [status,file_list_a]=system(command);
   file_list_b=split(file_list_a);
   file_list=squeeze(file_list_b);
   nfile=size(file_list);
%
% Constants
   Ru=8.316;
   Rd=286.9;
   eps=0.61;
   molec_wt_o3=.0480;
   molec_wt_no2=.0460;
   molec_wt_so2=.0641;
   AvogN=6.02214e23;
   msq2cmsq=1.e4;
   P_std=1013.25;
   grav=9.8;
   cone_fac=.715567;   
%
% Convert DU to moles/m^2
   du2molpm2=4.4615e-4;
%
% Convert DU to molecules/m^2
   du2molcpm2=2.6867e20;
%
   win_secs_beg=whh_mn*60.*60. + wmm_mn*60. + wss_mn;
   win_secs_end=whh_mx*60.*60. + wmm_mx*60. + wss_mx;
%
% Print input data
%   fprintf('obs window str %d %d %d %d %d %d \n',wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn)
%   fprintf('obs window end %d %d %d %d %d %d \n',wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx)
%
% Read model grid
   lon_mdl=ncread(strcat(path_mdl,'/',file_mdl),'XLONG');
   lat_mdl=ncread(strcat(path_mdl,'/',file_mdl),'XLAT');
   delx=ncreadatt(strcat(path_mdl,'/',file_mdl),'/','DX');  
   cen_lat=ncreadatt(strcat(path_mdl,'/',file_mdl),'/','CEN_LAT');  
   cen_lon=ncreadatt(strcat(path_mdl,'/',file_mdl),'/','CEN_LON');
   if(cen_lon<0)
      cen_lon=cen_lon+360.;
   end
   truelat1=ncreadatt(strcat(path_mdl,'/',file_mdl),'/','TRUELAT1');  
   truelat2=ncreadatt(strcat(path_mdl,'/',file_mdl),'/','TRUELAT2');  
   moad_cen_lat=ncreadatt(strcat(path_mdl,'/',file_mdl),'/','MOAD_CEN_LAT');  
   stand_lon=ncreadatt(strcat(path_mdl,'/',file_mdl),'/','STAND_LON');  
   pole_lat=ncreadatt(strcat(path_mdl,'/',file_mdl),'/','POLE_LAT');  
   pole_lon=ncreadatt(strcat(path_mdl,'/',file_mdl),'/','POLE_LON');
%
% Process satellite data
   for ifile=1:nfile
      file_in=char(file_list(ifile));
      fprintf('%d %s \n',ifile,file_in);
      if(isempty(file_in))
         continue
      end
      indx=strfind(file_in,file_pre)-1;
      if(isempty(indx))
         continue
      end
      file_hh=str2double(file_in(indx+30:indx+31));
      file_mm=str2double(file_in(indx+32:indx+33));
      file_ss=0;
      file_str_secs=file_hh*60.*60. + file_mm*60. + file_ss;
%       
      if(win_secs_end<file_str_secs)
%         fprintf('APM: if win end secs %d < file str secs %d, then skip file \n', ...
%         win_secs_end,file_str_secs);
         continue
      end
      fprintf('READ OMI DATA \n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read OMI data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% date data
      field='/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/';
      day=h5readatt(file_in,field,'GranuleDay');
      month=h5readatt(file_in,field,'GranuleMonth');
      year=h5readatt(file_in,field,'GranuleYear');
% amf(npxl,nscan) (60, 1644) (nodim)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/AirMassFactor';
      amf=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');
      amf(:,:)=amf(:,:)*scalef;
      dims=size(amf);
      npxl=dims(1);
      nscan=dims(2);
      clear dims
% climo_lev(npxl,nscan,nlay) (hPa)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ClimatologyLevels';
      climo_lev=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');
      climo_lev(:,:,:)=climo_lev(:,:,:)*scalef;
      dims=size(climo_lev);
      nlay=dims(3);
      nlev=nlay+1;
% hcho_col(npxl,nscan) (molec/cm^2)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ColumnAmount';
      hcho_col=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      hcho_col(:,:)=hcho_col(:,:)*scalef;
% hcho_col_err(npxl,nscan) (molec/cm^2)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ColumnUncertainty';
      hcho_col_err=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      hcho_col_err(:,:)=hcho_col_err(:,:)*scalef;
% scat_wt(npxl,nscan,nlay) (none) (60, 1644, 47) (nodim)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ScatteringWeights';
      scat_wt=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      scat_wt(:,:,:)=scat_wt(:,:,:)*scalef;
% lat(npxl,nscan)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/Latitude';
      lat=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      lat(:,:)=lat(:,:)*scalef;;  
% lon(npxl,nscan)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/Longitude';
      lon=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      lon(:,:)=lon(:,:)*scalef;;  
      for i=1:npxl
         for j=1:nscan
            if(lon(i,j)<0.)
      	       lon(i,j)=lon(i,j)+360.;
            end
         end
      end
% zenang(npxl,nscan) (deg)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/SolarZenithAngle';
      zenang=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      zenang(:,:)=zenang(:,:)*scalef;;  
% time(6,nscan)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/TimeUTC';
      time=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      title=h5readatt(file_in,field,'Title');
      units=h5readatt(file_in,field,'Units');
      time(:,:)=time(:,:)*scalef;
%
% xtrk_flg(npxl,nscan) (none)
      field='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/XtrackQualityFlags';
      xtrk_flg=h5read(file_in,field);
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');
      xtrk_flg(:,:)=xtrk_flg(:,:)*scalef;
%
% Loop through OMI data
      [windate_min]=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2010));
      [windate_max]=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2010));
      ocnt=0;
      icnt=0;
      for iscan=1:nscan
         yyyy_omi=single(time(1,iscan));
         mn_omi=single(time(2,iscan));
         dy_omi=single(time(3,iscan));
         hh_omi=single(time(4,iscan));
         mm_omi=single(time(5,iscan));
         ss_omi=single(time(6,iscan));
         omidate=single(convert_time_ref(yyyy_omi,mn_omi,dy_omi,hh_omi,mm_omi,ss_omi,2010));
%
% Check time
         if(omidate<windate_min | omidate>windate_max)
%            fprintf('min %d date %d max %d \n',windate_min,omidate,windate_max)
            continue
         end
         for ipxl=1:npxl
%
% QA/AC
            if((ipxl>=1 & ipxl<=5) | (ipxl>=56 & ipxl<=60))
%              fprintf('APM: ipxl issue \n') 
               continue
	    end
%
%	    if(bitand(vcd_flg(ipxl,iscan),1)~=0 | xtrk_flg(ipxl,iscan)~=0 | ...
%	    xtrk_flg(ipxl,iscan)~=255 | zenang(ipxl,iscan)>=80.)
%	       fprintf('APM: vcd or xtrk flag issue \n') 
%               continue
%	    end
%
%	    if(cld_rad_frac(ipxl,iscan)>=0.5 | terr_refl(ipxl,iscan)>=0.3)
%	       fprintf('APM: cld_frac or terr_refl issue \n') 
%               continue
%	    end
%
	    if(isnan(hcho_col(ipxl,iscan)) | hcho_col(ipxl,iscan)<=0)
%	       fprintf('APM: hcho_col  issue \n')
               continue
            end
	    if(isnan(hcho_col_err(ipxl,iscan)) | hcho_col_err(ipxl,iscan)<=0)
%	       fprintf('APM: hcho_col_err  issue \n')
               continue
            end

	    if(hcho_col(ipxl,iscan)>1.e18 | hcho_col_err(ipxl,iscan)>1.e18)
	       fprintf('APM: hcho is too large %d %d \n',hcho_col(ipxl,iscan),hcho_col_err(ipxl,iscan))
               continue
            end
%
%	    if(isnan(col_amt_trop(ipxl,iscan)) | col_amt_trop(ipxl,iscan)<=0)
%	       fprintf('APM: col_amt_trop issue issue \n') 
%               continue
%            end
%
%	    if(isnan(col_amt_trop_std(ipxl,iscan)) | col_amt_trop_std(ipxl,iscan)<=0)
%	       fprintf('APM: col_amt_trop_err issue \n') 
%               continue
%            end
%	    fprintf('APM: processing obs number %d %d \n ',iscan, ipxl)
%
% Check domain
% Input grid needs to be in degrees
% X coordinate is [0 to 360]
%		 
	    x_obser=lon(ipxl,iscan);
            y_obser=lat(ipxl,iscan);
            if(x_obser<0.)
	       x_obser=360.+x_obser;
            end
%
	    xmdl_sw=lon_mdl(1,1);
	    if(xmdl_sw<0.)
	       xmdl_sw=xmdl_sw+360.;
            end
%
            xmdl_mx=lon_mdl(nx_mdl,ny_mdl);
            if(xmdl_mx<0.)
               xmdl_mx=xmdl_mx+360.;
            end
%
% APM: Need to get this info from model
            [xi,xj]=w3fb13(y_obser,x_obser,lat_mdl(1,1), ...
	    xmdl_sw,delx,cen_lon,truelat1,truelat2);
            i_min = round(xi);
            j_min = round(xj);
            reject = 0;
%
% Check lower bounds
            if(i_min<1 & round(xi)==0)
	       i_min=1;
            elseif(i_min<1 & fix(xi)<0)
   	       i_min=-9999;
               j_min=-9999;
               reject=1;
            end
            if(j_min<1 & round(xj)==0)
               j_min=1;
            elseif (j_min<1 & fix(xj)<0)
               i_min=-9999;
               j_min=-9999;
               reject=1;
            end
%
% Check upper bounds
            if(i_min>nx_mdl & fix(xi)==nx_mdl)
               i_min=nx_mdl;
            elseif (i_min>nx_mdl & fix(xi)>nx_mdl)
               i_min=-9999;
               j_min=-9999;
               reject=1;
            end
            if(j_min>ny_mdl & fix(xj)==ny_mdl)
	       j_min=ny_mdl;
            elseif (j_min>ny_mdl & fix(xj)>ny_mdl)
               i_min=-9999;
               j_min=-9999;
               reject=1;
            end
            if(reject==1)
%               fprintf('x_mdl_min, x_obs, x_mdl_max: %6.2f %6.2f %6.2f \n',xmdl_sw, ...
%               x_obser,xmdl_mx)
%               fprintf('y_mdl_min, y_obs, y_mdl_max: %6.2f %6.2f %6.2f \n',lat_mdl(1,1), ...
%               y_obser,lat_mdl(nx_mdl,ny_mdl))
%               fprintf('i_min %d j_min %d \n',i_min,j_min)
               continue
            end
            if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
               fprintf('NO REJECT: i_min %d j_min %d \n',i_min,j_min)
               continue
            end
%
% Transfer pressures
            for ilay=1:nlay
               prs_lay(ilay)=climo_lev(ipxl,iscan,ilay);
            end
%
% Save data to ascii file
            icnt=icnt+1;
            fprintf(fid,'OMI_NO2_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_omi, ...
	    mn_omi,dy_omi,hh_omi,mm_omi,ss_omi);
	    fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,iscan),lon(ipxl,iscan));
            fprintf(fid,'%d %d \n',nlay,nlev);
            fprintf(fid,'%14.8g \n',amf(ipxl,iscan));
            fprintf(fid,'%14.8g %14.8g \n',hcho_col(ipxl,iscan), ...
            hcho_col_err(ipxl,iscan));
            fprintf(fid,'%14.8g \n',zenang(ipxl,iscan));
            fprintf(fid,'%14.8g ',scat_wt(ipxl,iscan,1:nlay));
            fprintf(fid,'\n');
 	    fprintf(fid,'%14.8g ',prs_lay(1:nlay));
            fprintf(fid,'\n');
         end
      end
      clear amfstrat amfstrat_clr amfstrat_cld amftrop amftrop_clr 
      clear cld_frac cld_prs cld_rad_frac col_amt col_amt_std 
      clear col_amt_trop col_amt_trop_std scat_wt scat_wt_prs 
      clear slnt_col_amt slnt_col_amt_std prs_trop vcd_flg 
      clear xtrk_flg lat lon zenang time terr_refl
   end
end

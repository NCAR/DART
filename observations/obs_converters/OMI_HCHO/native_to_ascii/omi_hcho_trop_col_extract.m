function omi_hcho_trop_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
   day_secs_beg=whh_mn*60.*60. + wmm_mn*60. + wss_mn;
   day_secs_end=whh_mx*60.*60. + wmm_mx*60. + wss_mx;
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
      if(isempty(file_in))
         continue
      end
      indx=strfind(file_in,file_pre)-1;
      if(isempty(indx))
         continue
      end
      file_mm=str2double(file_in(indx+31:indx+32));
      file_secs=file_mm*60.;
%       
      fprintf('%d %s \n',ifile,file_in);
      fprintf('file str secs %d cycle vend secs %d \n',file_secs,day_secs_end);
%      if(day_secs_end<file_secs)
%         continue
%      end
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
      field='/HDFEOS/SWATHS/ColumnAmountNO2/';
      zgrid=h5readatt(file_in,field,'VerticalCoordinate');
% amfstrat(pixel,scanline)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/AmfStrat';
      amfstrat=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      amfstrat(:,:)=amfstrat(:,:)*scalef;
% amfstrat_clr(pixel,scanline)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/AmfStratClear';
      amfstrat_clr=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      amfstrat_clr(:,:)=amfstrat_clr(:,:)*scalef;
% amfstrat_cld(pixel,scanline)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/AmfStratCloudy';
      amfstrat_cld=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      amfstrat_cld(:,:)=amfstrat_cld(:,:)*scalef;
% amftrop(pixel,scanline) (none)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/AmfTrop';
      amftrop=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      amftrop(:,:)=amftrop(:,:)*scalef;
% amftrop_clr(pixel,scanline) (none)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/AmfTropClear';
      amftrop_clr=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      amftrop_clr(:,:)=amftrop_clr(:,:)*scalef;
% amftrop_cld(pixel,scanline) (none)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/AmfTropCloudy';
      amftrop_cld=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      amftrop_cld(:,:)=amftrop_cld(:,:)*scalef;
% cld_frac(pixel,scanline)
      clear temp
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/CloudFraction';
      temp=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      cld_frac(:,:)=double(temp(:,:))*scalef;
% cld_prs(pixel,scanline)
      clear temp
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/CloudPressure';
      temp=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      cld_prs(:,:)=double(temp(:,:))*scalef;
% cld_rad_frac(pixel,scanline) (scalef=.001)
      clear temp
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/CloudRadianceFraction';
      temp=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      cld_rad_frac(:,:)=double(temp(:,:))*scalef;
% col_amt(pixel,scanline) (molec/cm^2)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2';
      col_amt=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      col_amt(:,:)=col_amt(:,:)*scalef;
% col_amt_std(pixel,scanline) (molec/cm^2)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Std';
      col_amt_std=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      col_amt_std(:,:)=col_amt_std(:,:)*scalef;
% col_amt_trop(pixel,scanline) (molec/cm^2)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Trop';
      col_amt_trop=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      col_amt_trop(:,:)=col_amt_trop(:,:)*scalef;
% col_amt_trop_std(pixel,scanline) (molec/cm^2)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2TropStd';
      col_amt_trop_std=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      col_amt_trop_std(:,:)=col_amt_trop_std(:,:)*scalef;
% scat_wt(layer,pixel,scanline) (none)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWeight';
      scat_wt=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      scat_wt(:,:,:)=scat_wt(:,:,:)*scalef;
      tmp=size(scat_wt);
      layer=tmp(1);
      pixel=tmp(2);
      scanline=tmp(3);
      level=layer+1;
% scat_wt_prs(layer,pixel,scanline) (hPa) (bottom to top)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWtPressure';
      scat_wt_prs=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      scat_wt_prs(:,:,:)=scat_wt_prs(:,:,:)*scalef;
% bottom to top
      prs_lev(level)=scat_wt_prs(layer) - (scat_wt_prs(layer-1) ...
      -scat_wt_prs(layer))/2.;
      for ilv=1:layer-1
         ilvv=layer-ilv+1;
         prs_lev(ilvv)=(scat_wt_prs(ilvv)+scat_wt_prs(ilvv-1))/2.;
      end
      prs_lev(1)=scat_wt_prs(1) + (scat_wt_prs(1)-scat_wt_prs(2))/2.;
      prs_lay(1:layer)=scat_wt_prs(1:layer);
% slnt_col_amt(pixel,scanline) (molec/cm^2)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/SlantColumnAmountNO2';
      slnt_col_amt=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      slnt_col_amt(:,:)=slnt_col_amt(:,:)*scalef;
% slnt_col_amt_std(pixel,scanline) (molec/cm^2)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/SlantColumnAmountNO2Std';
      slnt_col_amt_std=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      slnt_col_amt_std(:,:)=slnt_col_amt_std(:,:)*scalef;
% prs_trop(pixel,scanline)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/TropopausePressure';
      prs_trop=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      prs_trop(:,:)=prs_trop(:,:)*scalef;  
% vcd_flg(pixel,scanline) (none)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/VcdQualityFlags';
      vcd_flg=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      vcd_flg(:,:)=vcd_flg(:,:)*scalef;
% terr_refl(pixel,scanline) (none)
      clear temp
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/TerrainReflectivity';
      temp=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      terr_refl(:,:)=double(temp(:,:))*scalef;
% xtrk_flg(pixel,scanline) (none)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/XTrackQualityFlags';
      xtrk_flg=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      xtrk_flg(:,:)=xtrk_flg(:,:)*scalef;
% lat(pixel,scanline)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Latitude';
      lat=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      lat(:,:)=lat(:,:)*scalef;
% lon(pixel,scanline)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Longitude';
      lon=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      lon(:,:)=lon(:,:)*scalef;
      for i=1:pixel
         for j=1:scanline
            if(lon(i,j)<0.)
      	       lon(i,j)=lon(i,j)+360.;
            end
         end
      end
% zenang(pixel,scanline) (deg)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/SolarZenithAngle';
      zenang=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      units=h5readatt(file_in,field,'Units');
      zenang(:,:)=zenang(:,:)*scalef;
% time(scanline) TAI-93 time (secs)
      field='/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Time';
      time=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');
      title=h5readatt(file_in,field,'Title');
      units=h5readatt(file_in,field,'Units');
      time(:)=time(:)*scalef;
      time(:)=time(:)-37;
%
% Loop through OMI data
      clear temp
      [temp,rc]=time_tai93(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn);
      windate_min=single(temp);
      clear temp
      [temp,rc]=time_tai93(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx);
      windate_max=single(temp);
      clear temp
      [temp,rc]=time_tai93(year,month,day,0,0,0);
      base_omi=single(temp);
%      fprintf('%d \n',single(time(1)))
%      fprintf('%d %d %d \n',windate_min,base_omi,windate_max)
%   
      yyyy_omi=double(year);
      mn_omi=double(month);
      dy_omi=double(day);
      ocnt=0;
      icnt=0;
      for ilin=1:scanline
         secs_day(ilin)=time(ilin)-base_omi;
         hh_omi=double(idivide(int32(secs_day(ilin)),3600));
         mm_omi=double(idivide(mod(int32(secs_day(ilin)),3600),60));
         ss_omi=double(int32(secs_day(ilin))-int32(hh_omi*3600+mm_omi*60));
         omidate=single(time(ilin));
         if(int32(hh_omi)>23 | int32(mm_omi)>59 | int32(ss_omi)>59)
            fprintf('yr %d mn %d dy %d hh %d mm %d ss %d \n', ...
	    year,month,day,int32(hh_omi),int32(mm_omi),int32(ss_omi))
            [yr_nw,mn_nw,dy_nw,hh_nw,mm_nw,ss_nw,rc]=incr_time(year, ...
            month,day,int32(hh_omi),int32(mm_omi),int32(ss_omi));
            year=yr_nw;
            month=mn_nw;
            day=dy_nw;
            hh_omi=hh_nw;
            mm_omi=mm_nw;
            ss_omi=ss_nw;
%            fprintf('yr %d mn %d dy %d hh %d mm %d ss %d \n \n', ...
%   	    year,month,day,hh_omi,mm_omi,ss_omi)
	 end 
%
% Check time
%         fprintf('min %d date %d max %d \n',windate_min,omidate,windate_max)
         if(omidate<windate_min | omidate>windate_max)
            continue
         end
         for ipxl=1:pixel
%
% QA/AC
            if((ipxl>=1 & ipxl<=5) | (ipxl>=56 & ipxl<=60))
	       fprintf('APM: ipxl issue \n') 
	       continue
	    end
%
%	    if(bitand(vcd_flg(ipxl,ilin),1)~=0 | xtrk_flg(ipxl,ilin)~=0 | ...
%	    xtrk_flg(ipxl,ilin)~=255 | zenang(ipxl,ilin)>=80.)
%	       fprintf('APM: vcd or xtrk flag issue \n') 
%               continue
%	    end
%
	    if(cld_rad_frac(ipxl,ilin)>=0.5 | terr_refl(ipxl,ilin)>=0.3)
	       fprintf('APM: cld_frac or terr_refl issue \n ') 
               continue
	    end
%
	    if(isnan(slnt_col_amt(ipxl,ilin)) | slnt_col_amt(ipxl,ilin)<=0)
	       fprintf('APM: slnt_col_amt  issue \n ') 
               continue
            end
%
	    if(isnan(col_amt_trop(ipxl,ilin)) | col_amt_trop(ipxl,ilin)<=0)
	       fprintf('APM: col_amt_trop issue issue \n ') 
               continue
            end
%
	    if(isnan(col_amt_trop_std(ipxl,ilin)) | col_amt_trop_std(ipxl,ilin)<=0)
	       fprintf('APM: col_amt_trop_err issue \n ') 
               continue
            end
%	    fprintf('APM: processing obs number %d %d \n ',ilin, ipxl)
%
% Check domain
% Input grid needs to be in degrees
% X coordinate is [0 to 360]
%		 
	    x_obser=lon(ipxl,ilin);
            y_obser=lat(ipxl,ilin);
            if(x_obser<0.)
	       x_obser=360.+x_obser;
            end
%
	    xmdl_sw=lon_mdl(1,1);
	    if(xmdl_sw<0.)
	       xmdl_sw=xmdl_sw+360.;
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
	       continue
	    end
	    if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
	       continue
	    end
%
% Save data to ascii file
            icnt=icnt+1;
            fprintf(fid,'OMI_NO2_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_omi, ...
	    mn_omi,dy_omi,hh_omi,mm_omi,ss_omi);
	    fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,ilin),lon(ipxl,ilin));
            fprintf(fid,'%d %d \n',layer,level);
            fprintf(fid,'%14.8g %14.8g %14.8g \n',amfstrat(ipxl,ilin), ...
            amfstrat_clr(ipxl,ilin),amfstrat_cld(ipxl,ilin));
            fprintf(fid,'%14.8g %14.8g %14.8g \n',amftrop(ipxl,ilin), ...
            amftrop_clr(ipxl,ilin),amftrop_cld(ipxl,ilin));
            fprintf(fid,'%14.8g %14.8g %14.8g \n',cld_frac(ipxl,ilin), ...
            cld_prs(ipxl,ilin),cld_rad_frac(ipxl,ilin));
            fprintf(fid,'%14.8g %14.8g \n',col_amt(ipxl,ilin), ...
	    col_amt_std(ipxl,ilin));
            fprintf(fid,'%14.8g %14.8g \n',col_amt_trop(ipxl,ilin), ...
	    col_amt_trop_std(ipxl,ilin));
            fprintf(fid,'%14.8g %14.8g \n',slnt_col_amt(ipxl,ilin), ...
	    slnt_col_amt_std(ipxl,ilin));
            fprintf(fid,'%14.8g %14.8g \n',prs_trop(ipxl,ilin), ...
	    zenang(ipxl,ilin));
            fprintf(fid,'%14.8g ',scat_wt(1:layer,ipxl,ilin));
            fprintf(fid,'\n');
 	    fprintf(fid,'%14.8g ',prs_lev(1:level));
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

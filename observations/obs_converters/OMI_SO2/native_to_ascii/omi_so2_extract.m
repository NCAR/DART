function omi_so2_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
      file_mm=str2double(file_in(indx+29:indx+32));
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
% flg_snoice(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/AlgorithmFlag_SnowIce';
      flg_snoice=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% flg_rowanom(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/Flag_RowAnomaly';
      flg_rowanom=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% flg_saa(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/Flag_SAA';
      flg_saa=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% cld_frac(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/CloudFraction';
      cld_frac=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% cld_prs(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/CloudPressure';
      cld_prs=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% cld_rad_frac(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/CloudRadianceFraction';
      cld_rad_frac=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% col_amt(pixel,scanline) (Dobson Units)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/ColumnAmountSO2';
      col_amt=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      col_amt(:,:)=col_amt(:,:)*du2molcpm2/msq2cmsq;
% col_amt_pbl(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/ColumnAmountSO2_PBL';
      col_amt_pbl=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');
      range=h5readatt(file_in,field,'ValidRange');
      col_amt_pbl(:,:)=col_amt_pbl(:,:)*du2molcpm2/msq2cmsq;
% col_amt_stl(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/ColumnAmountSO2_STL';
      col_amt_stl=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% col_amt_trl(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/ColumnAmountSO2_TRL';
      col_amt_trl=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      col_amt_trl(:,:)=col_amt_trl(:,:)*du2molcpm2/msq2cmsq;
% col_amt_trm(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/ColumnAmountSO2_TRM';
      col_amt_trm=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      col_amt_trm(:,:)=col_amt_trm(:,:)*du2molcpm2/msq2cmsq;
% col_amt_tru(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/ColumnAmountSO2_TRU';
      col_amt_tru=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      col_amt_tru(:,:)=col_amt_tru(:,:)*du2molcpm2/msq2cmsq;
% layer_wt(layer,pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/GEOS5LayerWeight';
      layer_wt=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      tmp=size(layer_wt);
      layer=tmp(1);
      pixel=tmp(2);
      scanline=tmp(3);
      level=layer+1;
% prs_bot(layer) (hPa)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/LayerBottomPressure';
      prs_bot=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');
% layer_wt_pbl(layer,pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/PBLLayerWeight';
      layer_wt_pbl=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% rad_cld_frac(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/RadiativeCloudFraction';
      rad_cld_frac=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% scat_wt(layer,pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/ScatteringWeight';
      scat_wt=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');
% slnt_col_amt(pixel,scanline) (molecules/cm^2)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/SlantColumnAmountSO2';
      slnt_col_amt=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');
% calculate amf_total
      amf_total(:,:)=slnt_col_amt(:,:)./col_amt(:,:);
% lat(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Geolocation Fields/Latitude';
      lat=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% lon(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Geolocation Fields/Longitude';
      lon=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      for i=1:pixel
         for j=1:scanline
            if(lon(i,j)<0.)
      	       lon(i,j)=lon(i,j)+360.;
            end
         end
      end
% secs_day(scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Geolocation Fields/SecondsInDay';
      secs_day=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
% zenang(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Geolocation Fields/SolarZenithAngle';
      zenang=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      radcldfrc=h5read(file_in,field);
% time(scanline) TAI93 seconds
      field='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Geolocation Fields/Time';
      time=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');   
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      time(:)=time(:)-37;
%
% Define level pressures (top to bottom)
      prs_lev(1)=0.01;
      for ilv=2:level
         prs_lev(ilv)=prs_bot(ilv-1);
      end
%
% Define layer pressures
      for ilv=1:layer
         prs_lay(ilv)=(prs_lev(ilv)+prs_lev(ilv+1))/2.;
      end
%
% Loop through OMI data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2010));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2010));
      icnt=0;
      for ilin=1:scanline
         yyyy_omi=double(year);
         mn_omi=double(month);
         dy_omi=double(day);
         hh_omi=double(idivide(int32(secs_day(ilin)),3600));
         mm_omi=double(idivide(mod(int32(secs_day(ilin)),3600),60));
         ss_omi=double(int32(secs_day(ilin))-int32(hh_omi*3600+mm_omi*60));
         if(int32(hh_omi)>23 | int32(mm_omi)>59 | int32(ss_omi)>59)
            [yyyy_omi,mn_omi,dy_omi,hh_omi,mm_omi,ss_omi]=incr_time(yyyy_omi, ...
      	    mn_omi,dy_omi,hh_omi,mm_omi,ss_omi);
         end
         omidate=single(convert_time_ref(yyyy_omi,mn_omi,dy_omi,hh_omi,mm_omi,ss_omi,2010));
%
% Check time
         if(omidate<windate_min | omidate>windate_max)
            continue
         end
         omidate=single(convert_time_ref(yyyy_omi,mn_omi,dy_omi,hh_omi,mm_omi,ss_omi,2010));
%
% Check time
         if(omidate<windate_min | omidate>windate_max)
            continue
         end
%         for ipxl=1:pixel
         for ipxl=5:55
%
% QA/AC
%            fprintf(' \n');
%            fprintf('rowanom  %d saa %d crf %d \n',flg_rowanom(ipxl,ilin),flg_saa(ipxl,ilin),cld_rad_frac(ipxl,ilin));
%            fprintf('zenang  %d snoice %d amf %d \n',zenang(ipxl,ilin),flg_snoice(ipxl,ilin),amf_total(ipxl,ilin));
%            fprintf('lat_min  %d lat %d lat_max %d \n',lat_min,lat(ipxl,ilin),lat_max)
%            fprintf('lon_min  %d lon %d lon_max %d \n',lon_min,lon(ipxl,ilin),lon_max)

            if(flg_rowanom(ipxl,ilin)==1 | flg_saa(ipxl,ilin)==1 | ...
	       cld_rad_frac(ipxl,ilin)>=0.3 | zenang(ipxl,ilin)>=65.0 | ...
	       flg_snoice(ipxl,ilin)==2 | amf_total(ipxl,ilin)<=0.3)
               continue
            end
	    if(isnan(slnt_col_amt(ipxl,ilin)) | slnt_col_amt(ipxl,ilin)<=0)
               continue
            end
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
            fprintf(fid,'OMI_SO2_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_omi, ...
   	    mn_omi,dy_omi,hh_omi,mm_omi,ss_omi);
   	    fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,ilin),lon(ipxl,ilin));
            fprintf(fid,'%d %d \n',layer,level);
            fprintf(fid,'%14.8g %14.8g %14.8g \n',col_amt(ipxl,ilin), ...
            col_amt_pbl(ipxl,ilin),col_amt_stl(ipxl,ilin));
            fprintf(fid,'%14.8g %14.8g %14.8g \n',cld_frac(ipxl,ilin), ...
            cld_prs(ipxl,ilin),cld_rad_frac(ipxl,ilin));
            fprintf(fid,'%14.8g \n',rad_cld_frac(ipxl,ilin));
            fprintf(fid,'%14.8g \n',slnt_col_amt(ipxl,ilin));
            fprintf(fid,'%14.8g \n',zenang(ipxl,ilin));
            fprintf(fid,'%14.8g ',scat_wt(1:layer,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',prs_lev(1:level));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',layer_wt(1:layer,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',layer_wt_pbl(1:layer,ipxl,ilin));
            fprintf(fid,'\n');

%for k=1:layer
%   if(scat_wt(k,ipxl,ilin)==0)
%      fprintf('%d %d %d \n',ipxl,ilin,k)
%      fprintf('%14.8g',scat_wt(k:layer,ipxl,ilin))
%      fprintf(fid,'\n');
%      return	       
%   end
%end
         end
      end
      clear flg_snoice cld_frac cld_prs cld_rad_frac col_amt col_amt_pbl
      clear col_amt_stl col_amt_trl col_amt_trm col_amt_tru layer_wt
      clear prs_bot layer_wt_pbl rad_cld_frac scat_wt slnt_col_amt lat lon
      clear secs_day zenag time amf_total
      clear flg_rowanom flg_saa layer_wt_pbl
   end
end

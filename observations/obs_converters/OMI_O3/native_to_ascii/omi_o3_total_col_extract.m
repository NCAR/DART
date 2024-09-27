function omi_o3_total_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
   xcnt1=0;
   xcnt2=0;
   xcnt3=0;
   xcnt4=0;
   xcnt5=0;
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
      fprintf('min %d sec %d \n',file_mm,file_secs) 
%       
      fprintf('%d %s \n',ifile,file_in);
      fprintf('file str secs %d cycle end secs %d \n',file_secs,day_secs_end);
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
      field='/HDFEOS/SWATHS/OMI Column Amount O3/';
      ntimes=h5readatt(file_in,field,'NumTimes');
      npixel=h5readatt(file_in,field,'NumTimesSmallPixel');
      zgrid=h5readatt(file_in,field,'VerticalCoordinate');
% prior_lay(layer,pixel,scanline) (Dobson Units)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields/APrioriLayerO3';
      prior_lay=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      prior_lay(:,:,:)=prior_lay(:,:,:)*scalef;
      prior_lay(:,:,:)=prior_lay(:,:,:)*du2molpm2;
      tmp=size(prior_lay);
      layer=tmp(1);
      pixel=tmp(2);
      scanline=tmp(3);
      level=layer+1;
% cld_prs(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields/CloudPressure';
      cld_prs=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      cld_prs(:,:)=cld_prs(:,:)*scalef;
% col_amt(pixel,scanline) (Dobson Units)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields/ColumnAmountO3';
      col_amt=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');
      col_amt(:,:)=col_amt(:,:)*scalef;
      col_amt(:,:)=col_amt(:,:)*du2molpm2;
% qual_flg(pixel,scanline) (None)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields/QualityFlags';
      qual_flg=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');
      qual_flg(:,:)=qual_flg(:,:)*scalef;
% xtrk_flg(pixel,scanline) (None)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Geolocation Fields/XTrackQualityFlags';
      xtrk_flg=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');
      xtrk_flg(:,:)=xtrk_flg(:,:)*scalef;
% rad_cld_frc(pixel,scanline) (None)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields/RadiativeCloudFraction';
      rad_cld_frc=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      rad_cld_frc(:,:)=rad_cld_frc(:,:)*scalef;
% avgk_lay(layer,pixel,scanline) (None)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields/LayerEfficiency';
      avgk_lay=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      title=h5readatt(file_in,field,'Title');
      defn=h5readatt(file_in,field,'UniqueFieldDefinition');
      avgk_lay(:,:,:)=avgk_lay(:,:,:)*scalef;
% lat(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Geolocation Fields/Latitude';
      lat=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      lat(:,:)=lat(:,:)*scalef;
% lon(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Geolocation Fields/Longitude';
      lon=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      lon(:,:)=lon(:,:)*scalef;
      for i=1:pixel
         for j=1:scanline
            if(lon(i,j)<0.)
      	       lon(i,j)=lon(i,j)+360.;
            end
         end
      end
% secs_day(scanline)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Geolocation Fields/SecondsInDay';
      secs_day=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      secs_day(:)=secs_day(:)*scalef;
% zenang(pixel,scanline) (deg)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Geolocation Fields/SolarZenithAngle';
      zenang=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      zenang(:,:)=zenang(:,:)*scalef;
% time(scanline)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Geolocation Fields/Time';
      time=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');  
      range=h5readatt(file_in,field,'ValidRange');  
      time(:)=time(:)*scalef;
      time(:)=time(:)-37;
% sfc_prs(pixel,scanline)
      field='/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields/TerrainPressure';
      sfc_prs=h5read(file_in,field);
      missing=h5readatt(file_in,field,'MissingValue');  
      offset=h5readatt(file_in,field,'Offset');  
      scalef=h5readatt(file_in,field,'ScaleFactor');  
      units=h5readatt(file_in,field,'Units');
      range=h5readatt(file_in,field,'ValidRange');  
      sfc_prs(:)=sfc_prs(:)*scalef;
%
% Define OMI vertical pressure grid (hPa) (bottom to top)
      del_lnpr=log(2.0);
      lnpr(1)=0.;
      prs_lev(:,:,1)=exp(lnpr(1))*sfc_prs(:,:);
      for k=2:level
%         if(k==level)
%      	     prs_lev(:,:,k)=2.e-10*sfc_prs(:,:);
%            continue
%         end
         lnpr(k)=lnpr(k-1)-del_lnpr;
         prs_lev(:,:,k)=exp(lnpr(k))*sfc_prs(:,:);
      end
      for k=1:layer
         prs_lay(:,:,k)=(prs_lev(:,:,k)+prs_lev(:,:,k+1))/2.;
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
%	    fprintf('%d %d %d %d %d %d  \n',wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn)
%	    fprintf('%d %d %d %d %d %d  \n',yyyy_omi,mn_omi,dy_omi,hh_omi,mm_omi,ss_omi)
%	    fprintf('%d %d %d %d %d %d  \n',wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx)
%	    fprintf('\n')
% 	    fprintf('%d %d %d \n',windate_min,omidate,windate_max)
	    continue
         end
         for ipxl=1:pixel
%
% QA/AC	(qual_flg ~11000; xtrl_flg ~0;
	    if(zenang(ipxl,ilin) >= 85. | qual_flg(ipxl,ilin)~=0 | ...
	       xtrk_flg(ipxl,ilin)~=0)
 	       xcnt1=xcnt1+1;
	       continue
            end
	    if(ipxl==42 | ipxl==43 | ipxl==44 | ...
	      ipxl==45 | ipxl==46)
               continue
            end
	    if(isnan(col_amt(ipxl,ilin)) | col_amt(ipxl,ilin)<=0)
	      xcnt2=xcnt2+1;
               continue
            end
            if(isnan(prior_lay(layer,ipxl,ilin)) | prior_lay(layer,ipxl,ilin)<=0.)
	      xcnt3=xcnt3+1;
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
%            if(i_min<1 & round(xi)==0)
%	       i_min=1;
%            elseif(i_min<1 & fix(xi)<0)
%   	       i_min=-9999;
%               j_min=-9999;
%               reject=1;
%            end
%            if(j_min<1 & round(xj)==0)
%               j_min=1;
%            elseif (j_min<1 & fix(xj)<0)
%               i_min=-9999;
%               j_min=-9999;
%               reject=1;
%            end
%
% Check upper bounds
%            if(i_min>nx_mdl & fix(xi)==nx_mdl)
%               i_min=nx_mdl;
%            elseif (i_min>nx_mdl & fix(xi)>nx_mdl)
%               i_min=-9999;
%               j_min=-9999;
%               reject=1;
%            end
%            if(j_min>ny_mdl & fix(xj)==ny_mdl)
%	       j_min=ny_mdl;
%            elseif (j_min>ny_mdl & fix(xj)>ny_mdl)
%               i_min=-9999;
%               j_min=-9999;
%               reject=1;
%            end
            if(reject==1)
	      xcnt4=xcnt4+1;	      
	       continue
	    end
	    if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
	      xcnt5=xcnt5+1;
	       continue
	    end
%
% Save data to ascii file
            icnt=icnt+1;
            fprintf(fid,'OMI_TOMS_O3_Obs: %d %d %d  \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_omi, ...
	    mn_omi,dy_omi,hh_omi,mm_omi,ss_omi);
	    fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,ilin),lon(ipxl,ilin));
            fprintf(fid,'%d %d \n',layer,level);
   	    fprintf(fid,'%14.8g ',prs_lev(ipxl,ilin,1:level));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',avgk_lay(1:layer,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',prior_lay(1:layer,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g \n',col_amt(ipxl,ilin));
         end
      end
      clear prior_lay cld_prs col_amt rad_cld_frac avgk_lay 
      clear lat lon secs_day zenang time prs_lev prs_lay
   end
   fprintf('%d %d %d %d %d \n',xcnt1,xcnt2,xcnt3,xcnt4,xcnt5)
end
%

function tropomi_o3_profile_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
% Get file list and number of files
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
   molec_wt_co=.0480;
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
   day_secs_beg=whh_mn*60.*60. + wmm_mn*60. + wss_mn;
   day_secs_end=whh_mx*60.*60. + wmm_mx*60. + wss_mx;
%
% Print input data
   fprintf('obs window str %d %d %d %d %d %d \n',wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn)
   fprintf('obs window end %d %d %d %d %d %d \n',wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx)
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
      file_str_yy=str2double(file_in(indx+21:indx+24));
      file_str_mm=str2double(file_in(indx+25:indx+26));
      file_str_dd=str2double(file_in(indx+27:indx+28));
      file_str_hh=str2double(file_in(indx+30:indx+31));
      file_str_mn=str2double(file_in(indx+32:indx+33));
      file_str_ss=str2double(file_in(indx+34:indx+35));
      file_end_yy=str2double(file_in(indx+37:indx+40));
      file_end_mm=str2double(file_in(indx+41:indx+42));
      file_end_dd=str2double(file_in(indx+43:indx+44));
      file_end_hh=str2double(file_in(indx+46:indx+47));
      file_end_mn=str2double(file_in(indx+48:indx+49));
      file_end_ss=str2double(file_in(indx+50:indx+51));
      file_str_secs=file_str_hh*60.*60. + file_str_mn*60. + file_str_ss;
      file_end_secs=file_end_hh*60.*60. + file_end_mn*60. + file_end_ss;
      fprintf('%d %s \n',ifile,file_in);
      fprintf('file str %d cycle end %d \n',file_str_secs,day_secs_end);
      fprintf('file end %d cycle str %d \n',file_end_secs,day_secs_beg);
%
%%      if(file_str_secs>day_secs_end | file_end_secs<day_secs_beg)
%%         continue
%%      end
      fprintf('READ TROPOMI DATA \n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read TROPOMI data (profile data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read TROPOMI O3 data
%
% scanline
      field='/PRODUCT/scanline';   % 835
      temp=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');  
      scanline=size(temp);
      nscan=scanline(1);
% pixel
      field='/PRODUCT/ground_pixel';   % 77
      temp=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');  
      pixel=size(temp);
      npxl=pixel(1);
% time
      field='/PRODUCT/time';  % 1
      temp=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');  
      time=double(temp(1));
% nlay
      field='/PRODUCT/level';  % 33
      temp=ncread(file_in,field);
      level=size(temp);
      nlay=level(1);
      nlev=nlay+1;
% lat(npxl,nscan)
      field='/PRODUCT/latitude';
      lat=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% lon(npxl,nscan)
      field='/PRODUCT/longitude';
      lon=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      for ipxl=1:npxl
         for ilin=1:nscan
            if(lon(ipxl,ilin)<0)
      	       lon(ipxl,ilin)=lon(ipxl,ilin)+360.;
            end
         end
      end
% time_delta(nscan)
      field='/PRODUCT/delta_time';
      temp=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
      time_delta=double(temp(:,:))*1.e-3;
% time_utc(nscan) (APM: this fileld is blank)
      field='/PRODUCT/time_utc';
      time_utc=h5read(file_in,field);
% qa_value(npxl,nscan)
      field='/PRODUCT/qa_value';
      qa_value=ncread(file_in,field); 
      units=ncreadatt(file_in,field,'units');  
      scalef=ncreadatt(file_in,field,'scale_factor');
      offset=ncreadatt(file_in,field,'add_offset');
% o3_prof(nlay,npxl,nscan) (mol m-3)
      field='/PRODUCT/ozone_profile';
      o3_prof=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% o3_prof_err(nlay,npxl,nscan) (mol m-3)
      field='/PRODUCT/ozone_profile_precision';
      o3_prof_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% o3_prof_smth_err(nlay,npxl,nscan) (mol m-3)
      field='/PRODUCT/ozone_profile_smoothing';
      o3_prof_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% o3_total_col(npxl,nscan) (mol m-3)
      field='/PRODUCT/ozone_total_column';
      o3_total_col=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% o3_total_col_err(npxl,nscan) (mol m-3)
      field='/PRODUCT/ozone_total_column_precision';
      o3_total_col_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% o3_trop_col(npxl,nscan) (mol m-3)
      field='/PRODUCT/ozone_tropospheric_column';
      o3_trop_col=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% o3_trop_col_err(npxl,nscan) (mol m-3)
      field='/PRODUCT/ozone_tropospheric_column_precision';
      o3_trop_col_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% prs_lay(nlay,npxl,nscan) (Pa) (bottom to top)
      field='/PRODUCT/pressure';
      prs_lay=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      prs_lay=prs_lay/100.;
% solar_zenith_angle(npxl,nscan) (degrees)
      field='/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle';
      zenang=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% cov_err(nlay,nlay,npxl,nscan) (m)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/ozone_profile_error_covariance_matrix';
      cov_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% avgk_lay(nlay,nlay,npxl,nscan) (m)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/averaging_kernel';
      avgk_lay=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% dofs(npxl,nscan) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/degrees_of_freedom';
      dofs=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% prior_lay(nlay,npxl,nscan) (mol m-3)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/ozone_profile_apriori';
      prior_lay=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% prs_sfc(npxl,nscan)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure';
      prs_sfc=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% prs_trop(npxl,nscan)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/pressure_at_tropopause';
      prs_trop=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      fprintf('BEGIN DATA PROCESSING  \n')
%
% Loop through TROPOMI data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2010));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2010));
      icnt=0;
      for ilin=1:nscan
%
% Check time
      	 tropomidate=single(time+time_delta(ilin));
         [yyyy_tropomi,mn_tropomi,dy_tropomi,hh_tropomi,mm_tropomi, ...
         ss_tropomi]=invert_time_ref(tropomidate,2010);
         if(hh_tropomi<0 | mm_tropomi<0 | ss_tropomi<0)
            fprintf('%d %d %d \n',single(time),int64(time_delta(ipxl,ilin)),tropomidate)
            fprintf('%d %d %d %d %d %d \n',yyyy_tropomi,mn_tropomi,dy_tropomi, ...
            hh_tropomi,mm_tropomi,ss_tropomi);
            exit
         end
%
         if(tropomidate<windate_min | tropomidate>windate_max)
            continue
         end
         for ipxl=1:npxl
%
% QA/AC
%
            if(any(isnan(o3_prof(:,ipxl,ilin))) | any(o3_prof(:,ipxl,ilin)<=0))
               fprintf('FAILED QA/QC TEST \n')
               continue
            end
%
	    if(qa_value(ipxl,ilin)<0.50 | zenang(ipxl,ilin)>=80.0)
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
%              fprintf('FAILED DOMAIN TEST \n')
	       continue
	    end
	    if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
%              fprintf('FAILED DOMAIN TEST \n')
	       continue
	    end
%            fprintf('PASSED DOMAIN TEST \n')
%
% Save data to ascii file
            icnt=icnt+1;
            fprintf(fid,'TROPOMI_O3_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_tropomi, ...
            mn_tropomi,dy_tropomi,hh_tropomi,mm_tropomi,ss_tropomi);
            fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,ilin),lon(ipxl,ilin));
            fprintf(fid,'%d %d \n',nlay,nlev);
            fprintf(fid,'%14.8g ',prs_lay(1:nlay,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',avgk_lay(1:nlay,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',prior_lay(1:nlay,ipxl,ilin));
         end
      end
   end
end


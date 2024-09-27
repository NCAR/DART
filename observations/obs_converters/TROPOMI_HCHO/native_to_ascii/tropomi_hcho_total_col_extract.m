function tropomi_hcho_total_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
      time_start=ncreadatt(file_in,'/','time_coverage_start');
      time_end=ncreadatt(file_in,'/','time_coverage_end');  
      file_str_yy=str2double(time_start(1:4));
      file_str_mm=str2double(time_start(6:7));
      file_str_dd=str2double(time_start(9:10));
      file_str_hh=str2double(time_start(12:13));
      file_str_mn=str2double(time_start(15:16));
      file_str_ss=str2double(time_start(18:19));
      file_end_yy=str2double(time_end(1:4));
      file_end_mm=str2double(time_end(6:7));
      file_end_dd=str2double(time_end(9:10));
      file_end_hh=str2double(time_end(12:13));
      file_end_mn=str2double(time_end(15:16));
      file_end_ss=str2double(time_end(18:19));
      file_str_secs=file_str_hh*60.*60. + file_str_mn*60. + file_str_ss;
      file_end_secs=file_end_hh*60.*60. + file_end_mn*60. + file_end_ss;
%      fprintf('%d %s \n',ifile,file_in);
%      fprintf('file str %d cycle end %d \n',file_str_secs,day_secs_end);
%      fprintf('file end %d cycle str %d \n',file_end_secs,day_secs_beg);
%       
      if(file_str_secs>day_secs_end | file_end_secs<day_secs_beg)
         continue
      end
      fprintf('READ TROPOMI DATA \n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read TROPOMI data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% avgk_lay(nlay,npxl,nscan) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/averaging_kernel';
      avgk_lay=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      units=ncreadatt(file_in,field,'units');
      temp=size(avgk_lay);
      ntim=1;  
      nlay=temp(1);    % 34
      nlev=nlay+1;     % 35
      npxl=temp(2);    % 450
      nscan=temp(3);   % 4173
% slnt_col_mat(ncol,npxl,nscan) (mole/m^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/fitted_slant_columns';
      slnt_col_mat=double(ncread(file_in,field));
      indx_meaning=ncreadatt(file_in,field,'index_meaning');   
      long_name=ncreadatt(file_in,field,'long_name');   
      du_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_DU');   
      molec2cm2_conv_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_molecules_percm2');
      units=ncreadatt(file_in,field,'units');
% slnt_col_err_mat(ncol,npxl,nscan) (mole/m^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/fitted_slant_columns_precision';
      slnt_col_err_mat=double(ncread(file_in,field));
      indx_meaning=ncreadatt(file_in,field,'index_meaning');   
      long_name=ncreadatt(file_in,field,'long_name');   
      du_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_DU');   
      molec2cm2_conv_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_molecules_percm2');
      units=ncreadatt(file_in,field,'units');
% hcho_amf_clear(npxl,nscan) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/formaldehyde_clear_air_mass_factor';
      hcho_amf_clear=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      units=ncreadatt(file_in,field,'units');
% hcho_amf_cldy(npxl,nscan) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/formaldehyde_cloudy_air_mass_factor';
      hcho_amf_cldy=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      units=ncreadatt(file_in,field,'units');
% hcho_prof_prior(nlay,npxl,nscan)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/formaldehyde_profile_apriori';
      hcho_prof_prior=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      units=ncreadatt(file_in,field,'units');
% hcho_slnt_col_corr(npxl,nscan) (mole/m^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/formaldehyde_slant_column_corrected';
      hcho_slnt_col_corr=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      du_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_DU');   
      molec2cm2_conv_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_molecules_percm2');
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% hcho_slnt_col_corr_tru(npxl,nscan) (mole/m^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/formaldehyde_slant_column_corrected_trueness';
      hcho_slnt_col_corr_tru=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      du_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_DU');   
      molec2cm2_conv_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_molecules_percm2');
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% hcho_slnt_col_corr_flg(npxl,nscan) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/formaldehyde_slant_column_correction_flag';
      hcho_slnt_col_corr_flg=double(ncread(file_in,field));
      flg_meaning=ncreadatt(file_in,field,'flag_meanings');   
      flg_values=ncreadatt(file_in,field,'flag_values');
      long_name=ncreadatt(file_in,field,'long_name');   
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% hcho_amf_trop(npxl,nscan) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/formaldehyde_tropospheric_air_mass_factor';
      hcho_amf_trop=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      units=ncreadatt(file_in,field,'units');
% hcho_trop_col_correction(npxl,nscan)  (mole/m^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/formaldehyde_tropospheric_vertical_column_correction';
      hcho_trop_col_correction=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      du_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_DU');   
      molec2cm2_conv_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_molecules_percm2');
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% hcho_trop_col_tru(npxl,nscan) (mole/m^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/formaldehyde_tropospheric_vertical_column_trueness';
      hcho_trop_col_tru=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      du_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_DU');   
      molec2cm2_conv_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_molecules_percm2');
      units=ncreadatt(file_in,field,'units');
% solar_zenith_angle(npxl,nscan) (degrees)
      field='/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle';
      zenang=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% prs_sfc(npxl,nscan) (Pa)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure';
      prs_sfc=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% tm5_a(nlay) (Pa)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/tm5_constant_a';
      tm5_a=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% tm5_b(nlay) (none)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/tm5_constant_b';
      tm5_b=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
%
% Define TROPOMI vertical pressure grid (hPa) (bottom to top)
      for ipxl=1:npxl
         for ilin=1:nscan
            for ilv=1:nlay
               prs_lay(ilv,ipxl,ilin)=tm5_a(ilv)+tm5_b(ilv)* ...
               prs_sfc(ipxl,ilin);
#              if(prs_lay(ilv,ipxl,ilin)<.1)
#                 prs_lay(ilv,ipxl,ilin)=.1;
#              end
	    end
	    prs_lev(1,ipxl,ilin)=prs_sfc(ipxl,ilin)
	    for ilv=2:nlev
               prs_lev(ilv,ipxl,ilin)=2.*prs_lay(ilv-1,ipxl,ilin)-prs_lev(ilv-1,ipxl,ilin)
            end
        end
      end
      prs_lev=prs_lev/100.;
      prs_lay=prs_lay/100.;
% tm5_trop_indx(npxl,nscan) (none)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/tm5_tropopause_layer_index';
      tm5_trop_indx=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% time_delta(nplx,nscan) (millisecs)
      field='/PRODUCT/delta_time';
      time_delt=ncread(file_in,field);
      long_name=ncreadatt(file_in,field,'long_name');  
      units=ncreadatt(file_in,field,'units');
      time_delt=time_delt/1000.;
% hcho_trop_col(npxl,nscan) (mole/m^2)
      field='/PRODUCT/formaldehyde_tropospheric_vertical_column';
      hcho_trop_col=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      du_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_DU');   
      molec2cm2_conv_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_molecules_percm2');
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% hcho_trop_col_err(npxl,nscan) (mole/m^2)
      field='/PRODUCT/formaldehyde_tropospheric_vertical_column_precision';
      hcho_trop_col_err=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      du_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_DU');   
      molec2cm2_conv_conv=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_molecules_percm2');
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% lat(npxl,nscan)
      field='/PRODUCT/latitude';
      lat=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% layer(nlay)
      field='/PRODUCT/layer';
      layer=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      units=ncreadatt(file_in,field,'units');
% lon(npxl,nscan)
      field='/PRODUCT/longitude';
      lon=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
      for ipxl=1:npxl
         for ilin=1:nscan
            if(lon(ipxl,ilin)<0)
      	       lon(ipxl,ilin)=lon(ipxl,ilin)+360.;
            end
         end
      end
% qa_value(npxl,nscan)
      field='/PRODUCT/qa_value';
      qa_value=ncread(file_in,field); 
      offset=ncreadatt(file_in,field,'add_offset');
      long_name=ncreadatt(file_in,field,'long_name');   
      scalef=ncreadatt(file_in,field,'scale_factor');
      units=ncreadatt(file_in,field,'units');  
% scanline(nscan)
      field='/PRODUCT/scanline';
      scanline=ncread(file_in,field); 
      long_name=ncreadatt(file_in,field,'long_name');   
      units=ncreadatt(file_in,field,'units');  
% time(ntim) (seconds since 2010-01-01 00:00:00)
      field='/PRODUCT/time';
      time=ncread(file_in,field); 
      long_name=ncreadatt(file_in,field,'long_name');   
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% time_utc(nscan) (NETCDF and character data problem)
%      field='/PRODUCT/time_utc';
%      time_utc=ncread(file_in,field)
%      time_utc_cell=h5read(file_in,field);
%      time_utc=cell2mat(time_utc_cell);
%      time_utc(1)
%      long_name=ncreadatt(file_in,field,'long_name');   
%
      fprintf('BEGIN DATA PROCESSING \n')
%
% Loop through TROPOMI data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2010));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2010));
      icnt=0;
      for ipxl=1:npxl
         for ilin=1:nscan
            tropomidate=single(convert_time_ref(file_str_yy,file_str_mm,file_str_dd, ...
            0,0,0,2010))+time_delt(ipxl,ilin);
            [yyyy_tropomi,mn_tropomi,dd_tropomi,hh_tropomi,mm_tropomi,ss_tropomi]= ...
	    invert_time_ref(tropomidate,2010);
%            fprintf('windate_min %d \n',windate_min)
%            fprintf('tropomi_dat %d \n',tropomidate)
%            fprintf('windate_max %d \n',windate_max)
%
% Check time
            if(tropomidate<windate_min | tropomidate>windate_max)
               continue
            end
%
% QA/AC
%	    if(qa_value(ipxl,ilin)<0.50 | zenang(ipxl,ilin)>=80.0 | cld_rad_frac(ipxl,ilin) >=.5)
%               continue
%	    end
            if(isnan(hcho_trop_col(ipxl,ilin)) | hcho_trop_col(ipxl,ilin)<=0)
               continue
            end
%
            if(isnan(hcho_trop_col_err(ipxl,ilin)) | hcho_trop_col_err(ipxl,ilin)<=0)
               continue
            end
%
            reject=0;
            for ilay=1:nlay
               if (isnan(avgk_lay(ilay,ipxl,ilin)))
                  reject=1
                  break
               end
            end
            if(reject==1)
               continue
            end  
%
            reject=0;
            for ilay=1:nlay
               if (isnan(hcho_prof_prior(ilay,ipxl,ilin)))
                  reject=1
                  break
               end
            end
            if(reject==1)
               continue
            end  
%           fprintf('PASSED QA/QC TEST \n')
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
            elseif(i_min<1 & fix(xi)<=0)
   	       i_min=-9999;
               j_min=-9999;
               reject=1;
            end
            if(j_min<1 & round(xj)==0)
               j_min=1;
            elseif (j_min<1 & fix(xj)<=0)
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
% Save data to ascii file
	    icnt=icnt+1;
            fprintf(fid,'TROPOMI_HCHO_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_tropomi, ...
            mn_tropomi,dd_tropomi,hh_tropomi,mm_tropomi,ss_tropomi);
            fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,ilin),lon(ipxl,ilin));
            fprintf(fid,'%d %d \n',nlay,nlev);
            fprintf(fid,'%14.8g ',prs_lev(1:nlev,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',avgk_lay(1:nlay,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',hcho_prof_prior(1:nlay,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g %14.8g \n',hcho_trop_col(ipxl,ilin), ...
            hcho_trop_col_err(ipxl,ilin));
            fprintf(fid,'%14.8g \n',hcho_amf_trop(ipxl,ilin));
            fprintf(fid,'%d \n',tm5_trop_indx(ipxl,ilin));         
         end
      end
   end
end

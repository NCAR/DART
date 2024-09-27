function tropomi_ch4_cpsr_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
   molec_wt_no2=.0480;
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
% scanline
      field='/PRODUCT/scanline';
      temp=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');  
      long_name=ncreadatt(file_in,field,'long_name');  
      scanline=max(temp);
% pixel
      field='/PRODUCT/ground_pixel';
      temp=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');  
      long_name=ncreadatt(file_in,field,'long_name');  
      pixel=max(temp);
% time
%      field='/PRODUCT/ground_pixel';
%      temp=ncread(file_in,field);
%      units=ncreadatt(file_in,field,'units');  
%      long_name=ncreadatt(file_in,field,'long_name');  
%      time=max(temp);
% layer
      field='/PRODUCT/layer';
      layer_hgt=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');  
      standard_name=ncreadatt(file_in,field,'standard_name');  
      long_name=ncreadatt(file_in,field,'long_name');   
      tmp=size(layer_hgt);
      layer=tmp(1);
      level=layer+1;
% lat(pixel,scanline)
      field='/PRODUCT/latitude';
      lat=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      units=ncreadatt(file_in,field,'units');
      standard_name=ncreadatt(file_in,field,'standard_name');
% lon(pixel,scanline)
      field='/PRODUCT/longitude';
      lon=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      units=ncreadatt(file_in,field,'units');
      standard_name=ncreadatt(file_in,field,'standard_name');
      for ipxl=1:pixel
         for ilin=1:scanline
            if(lon(ipxl,ilin)<0)
      	       lon(ipxl,ilin)=lon(ipxl,ilin)+360.;
            end
         end
      end
% time_delta(scanline)
%      field='/PRODUCT/delta_time';
%      temp=ncread(file_in,field);
%      long_name=ncreadatt(file_in,field,'long_name');  
%      units=ncreadatt(file_in,field,'units');  
%      time_dela=double(temp)*1.e3;
% time_utc(scanline)
      field='/PRODUCT/time_utc';
      time_utc=h5read(file_in,field);
% qa_value(pixel,scanline)
      field='/PRODUCT/qa_value';
      qa_value=ncread(file_in,field); 
      units=ncreadatt(file_in,field,'units');  
      scalef=ncreadatt(file_in,field,'scale_factor');
      offset=ncreadatt(file_in,field,'add_offset');
      long_name=ncreadatt(file_in,field,'long_name');   
% col_amt_trop(pixel,scanline) (mol m-2)
      field='/PRODUCT/nitrogendioxide_tropospheric_column';
      col_amt_trop=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      standard_name=ncreadatt(file_in,field,'standard_name');
      long_name=ncreadatt(file_in,field,'long_name');   
% col_amt_trop_err(pixel,scanline) (mol m-2)
      field='/PRODUCT/nitrogendioxide_tropospheric_column_precision';
      col_amt_trop_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      standard_name=ncreadatt(file_in,field,'standard_name');
      long_name=ncreadatt(file_in,field,'long_name');   
% avgk_lay(layer,pixel,scanline) (none)
      field='/PRODUCT/averaging_kernel';
      avgk_lay=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% air_mass_factor_troposphere(pixel,scanline) (none)
      field='/PRODUCT/air_mass_factor_troposphere';
      amf_trop=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% air_mass_factor_total(pixel,scanline) (none)
      field='/PRODUCT/air_mass_factor_total';
      amf_total=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% indx_trop(pixel,scanline) (none)
      field='/PRODUCT/tm5_tropopause_layer_index';
      indx_trop=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% solar_zenith_angle(pixel,scanline) (degrees)
      field='/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle';
      zenang=double(ncread(file_in,field));
      long_name=ncreadatt(file_in,field,'long_name');   
      standard_name=ncreadatt(file_in,field,'standard_name');
      units=ncreadatt(file_in,field,'units');
% col_amt_strat(pixel,scanline) (mol m-2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_stratospheric_column';
      col_amt_strat=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
      standard_name=ncreadatt(file_in,field,'standard_name');
% col_amt_strat_err(pixel,scanline) (mol m-2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_stratospheric_column_precision';
      col_amt_strat_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
      standard_name=ncreadatt(file_in,field,'standard_name');
% col_amt_total(pixel,scanline) (mol m-2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_total_column';
      col_amt_total=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% col_amt_total_err(pixel,scanline) (mol m-2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_total_column_precision';
      col_amt_total_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% col_amt_summed_total(pixel,scanline) (mol m-2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_summed_total_column';
      col_amt_summed_total=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% col_amt_summed_total_err(pixel,scanline) (mol m-2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_summed_total_column_precision';
      col_amt_summed_total_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% slnt_col_amt(pixel,scanline) (mol m-2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_slant_column_density';
      slnt_col_amt=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% slnt_col_amt_err(pixel,scanline) (mol m-2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_slant_column_density_precision';
      slnt_col_amt_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% cld_rad_frac(pixel,scanline) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_radiance_fraction_nitrogendioxide_window';
      cld_rad_frac=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% dofs(pixel,scanline) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/degrees_of_freedom';
      dofs=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% amf_strat(pixel,scanline) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/air_mass_factor_stratosphere';
      amf_strat=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% amf_cld(pixel,scanline) (none)
%      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/air_mass_factor_cloudy';
%      amf_cld=double(ncread(file_in,field));
%      units=ncreadatt(file_in,field,'units');
%      long_name=ncreadatt(file_in,field,'long_name');   
% amf_clr(pixel,scanline) (none)
%      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/air_mass_factor_clear';
%      amf_clr=double(ncread(file_in,field));
%      units=ncreadatt(file_in,field,'units');
%      long_name=ncreadatt(file_in,field,'long_name');   
% tm5_a(vertices,layer)
      field='/PRODUCT/tm5_constant_a';
      tm5_a=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% tm5_b(vertices,layer)
      field='/PRODUCT/tm5_constant_b';
      tm5_b=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      long_name=ncreadatt(file_in,field,'long_name');   
% prs_sfc(pixels,scanline) (Pa)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure';
      prs_sfc=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      standard_name=ncreadatt(file_in,field,'standard_name');
      long_name=ncreadatt(file_in,field,'long_name');   
      fprintf('BEGIN DATA PROCESSING \n')
%
% Define TROPOMI vertical pressure grid (hPa) (bottom to top)
         for ipxl=1:pixel
            for ilin=1:scanline
               for ilv=1:layer
                  prs_lev(ilv,ipxl,ilin)=tm5_a(1,ilv)+tm5_b(1,ilv)* ...
                  prs_sfc(ipxl,ilin);
#                  if(prs_lev(ilv+1,ipxl,ilin)<.1)
#                     prs_lev(ilv+1,ipxl,ilin)=.1;
#                  end
                  prs_lay(ilv,ipxl,ilin)=(tm5_a(1,ilv)+tm5_b(1,ilv)* ...
                  prs_sfc(ipxl,ilin) + tm5_a(2,ilv)+tm5_b(2,ilv)* ...
                  prs_sfc(ipxl,ilin))/2.;
               end
               prs_lev(layer+1,ipxl,ilin)=tm5_a(2,layer)+tm5_b(2,layer)* ...
               prs_sfc(ipxl,ilin);
            end
         end
         prs_lev=prs_lev/100.;
         prs_lay=prs_lay/100.;
%
% Loop through TROPOMI data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2010));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2010));
      icnt=0;
      for ilin=1:scanline
         date_str=char(time_utc(ilin));
         yyyy_tropomi=str2double(extractBetween(date_str,1,4));
         mn_tropomi=str2double(extractBetween(date_str,6,7));
         dy_tropomi=str2double(extractBetween(date_str,9,10));
         hh_tropomi=str2double(extractBetween(date_str,12,13));
         mm_tropomi=str2double(extractBetween(date_str,15,16));
         ss_tropomi=str2double(extractBetween(date_str,18,26));
         if(int32(hh_tropomi)>23 | int32(mm_tropomi)>59 | ...
         int32(ss_tropomi)>59)
            [yyyy_tropomi,mn_tropomi,dy_tropomi,hh_tropomi, ...
            mm_tropomi,ss_tropomi]=incr_time(yyyy_tropomi, ...
      	 mn_tropomi,dy_tropomi,hh_tropomi,mm_tropomi,ss_tropomi);
         end
%         fprintf('obs date/time %d %d %d %d %d %d \n',yyyy_tropomi, ...
%         mn_tropomi,dy_tropomi,hh_tropomi,mm_tropomi,ss_tropomi)
         tropomidate=single(convert_time_ref(yyyy_tropomi,mn_tropomi, ...
         dy_tropomi,hh_tropomi,mm_tropomi,ss_tropomi,2010));
%         fprintf('windate_min %d \n',windate_min)
%         fprintf('tropomi_dat %d \n',tropomidate)
%         fprintf('windate_max %d \n',windate_max)
%
% Check time
         if(tropomidate<windate_min | tropomidate>windate_max)
            continue
         end
%         fprintf('PASSED DATE/TIME TEST \n')
         for ipxl=1:pixel
%
% QA/AC
	    if(qa_value(ipxl,ilin)<0.50 | zenang(ipxl,ilin)>=80.0 | cld_rad_frac(ipxl,ilin) >=.5)
               continue
	    end
            if(isnan(col_amt_trop(ipxl,ilin)) | col_amt_trop(ipxl,ilin)<=0)
               continue
            end
%            fprintf('PASSED QA/QC TEST \n')
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
%               fprintf('FAILED DOMAIN TEST \n')
	       continue
	    end
	    if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
%            fprintf('FAILED DOMAIN TEST \n')
	       continue
	    end
%            fprintf('PASSED DOMAIN TEST \n')
%
% Save data to ascii file
	    icnt=icnt+1;
            fprintf(fid,'TROPOMI_CH4_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_tropomi, ...
            mn_tropomi,dy_tropomi,hh_tropomi,mm_tropomi,ss_tropomi);
            fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,ilin),lon(ipxl,ilin));
            fprintf(fid,'%d %d \n',layer,level);
            fprintf(fid,'%14.8g ',prs_lev(1:level,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',avgk_lay(1:layer,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g %14.8g \n',col_amt_trop(ipxl,ilin), ...
            col_amt_trop_err(ipxl,ilin));
            fprintf(fid,'%14.8g \n',amf_trop(ipxl,ilin));
            fprintf(fid,'%d \n',indx_trop(ipxl,ilin));         
         end
      end
   end
end

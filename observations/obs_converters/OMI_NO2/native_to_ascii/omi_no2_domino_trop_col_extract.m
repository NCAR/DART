function omi_no2_domino_trop_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
   command=strcat('/usr/bin/rm'," ",'-rf'," ",fileout);
   [status]=system(command);
   fid=fopen(fileout,'w');
%
   command=strcat('/usr/bin/ls'," ",'-1'," ",filein,'*');
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
      time_start=ncreadatt(file_in,'/','time_coverage_start');
      time_ref=ncreadatt(file_in,'/','time_reference');
      time_end=ncreadatt(file_in,'/','time_coverage_end');
      ref_yy=str2double(time_ref(1:4));
      ref_mm=str2double(time_ref(6:7));
      ref_dd=str2double(time_ref(9:10));
      ref_hh=str2double(time_ref(12:13));
      ref_mn=str2double(time_ref(15:16));
      ref_ss=str2double(time_ref(18:19));
      ref_secs=single(convert_time_ref(ref_yy,ref_mm,ref_dd,ref_hh,ref_mn,ref_ss,1995));
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
      file_str_secs=file_str_hh*3600 + file_str_mn*60 + file_str_ss;
      file_end_secs=file_end_hh*3600 + file_end_mn*60 + file_end_ss;
      fprintf('%d %s \n',ifile,file_in);
      fprintf('If file_str_secs %d <= day_secs_end %d and \n',file_str_secs,day_secs_end);
      fprintf('   file_end_secs %d >= day_secs_beg %d then process data \n',file_end_secs,day_secs_beg);
%       
      if(file_str_secs>day_secs_end | file_end_secs<day_secs_beg)
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
      wfid=netcdf.open(file_in,'NC_NOWRITE');
      gfid=netcdf.inqNcid(wfid,'PRODUCT');
      dimid=netcdf.inqDimID(gfid,'scanline');
      [name,nscan]=netcdf.inqDim(gfid,dimid);              % 1644
      dimid=netcdf.inqDimID(gfid,'ground_pixel');
      [name,npxl]=netcdf.inqDim(gfid,dimid);               % 60
      dimid=netcdf.inqDimID(gfid,'corner');
      [name,ncnr]=netcdf.inqDim(gfid,dimid);               % 4
      dimid=netcdf.inqDimID(gfid,'time');
      [name,ntim]=netcdf.inqDim(gfid,dimid);               % 1
      dimid=netcdf.inqDimID(gfid,'polynomial_exponents');  
      [name,nply]=netcdf.inqDim(gfid,dimid);               % 5
      dimid=netcdf.inqDimID(gfid,'layer');
      [name,layer]=netcdf.inqDim(gfid,dimid);              % 34
      level=layer+1;
      dimid=netcdf.inqDimID(gfid,'nv');
      [name,nv]=netcdf.inqDim(gfid,dimid);                 % 2
%
% lat(npxl,nscan,ntim)
      field='/PRODUCT/latitude';
      lat=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% lon(npxl,nscan,ntim)
      field='/PRODUCT/longitude';
      lon=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
      for i=1:npxl
         for j=1:nscan
            if(lon(i,j)<0.)
      	       lon(i,j)=lon(i,j)+360.;
            end
         end
      end
%
% time(ntim)
      field='/PRODUCT/time';
      time=ncread(file_in,field);
% delta_time(nscan,ntim) (milliseconds)
      field='/PRODUCT/delta_time';
      delta_time=ncread(file_in,field);
      delta_time=delta_time/1000.;
      units=ncreadatt(file_in,field,'units');
% processing_error_flag(npxl,nscan,ntim) (8-bit unsigned integer; 0 - good, 1 - fail)
      field='/PRODUCT/processing_error_flag';
      proc_err_flag=ncread(file_in,field);
% snow_ice_flag(npxl,nscan,ntim) (8-bit unsigned integer; 0 - good, 1 - fail)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/snow_ice_flag';
      snow_ice_flag=ncread(file_in,field);
% no2_trop_col(npxl,nscan,ntim) (molec/cm^2)
      field='/PRODUCT/tropospheric_no2_vertical_column';
      no2_trop_col=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% no2_trop_col_err(npxl,nscan,ntim) (molec/cm^2)
      field='/PRODUCT/tropospheric_no2_vertical_column_uncertainty';
      no2_trop_col_err=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% avgk(layer,npxl,nscan,ntime) ( )
      field='/PRODUCT/averaging_kernel';
      avgk=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% amf_trop(npxl,nscan,ntime) ( )
      field='/PRODUCT/amf_trop';
      amf_trop=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% amf_total(npxl,nscan,ntime) ( )
      field='/PRODUCT/amf_total';
      amf_total=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% amf_geom(npxl,nscan,ntime) ( )
      field='/PRODUCT/amf_total';
      amf_geom=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% cloud_radiance_fraction_no2(npxl,nscan,ntime) ( )
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_radiance_fraction_no2';
      cld_rad_frac_no2=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% trop_indx(npxl,nscan,ntime) ( ) this is fill values only
      field='/PRODUCT/tm5_tropopause_layer_index';
      trop_indx=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% tm5_prs_a(nv,layer) (Pa) Convert to hPa (1 - lower interface, 2 - upper interface)
      field='/PRODUCT/tm5_pressure_level_a';
      tm5_prs_a=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
      tm5_prs_a=tm5_prs_a/100.;
% tm5_prs_b(nv,layer) ( )  (1 - lower interface, 2 - upper interface)
      field='/PRODUCT/tm5_pressure_level_b';
      tm5_prs_b=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% tm5_prs_sfc(npxl,nscan,ntime) (hPa)
      field='/PRODUCT/tm5_surface_pressure';
      tm5_prs_sfc=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
%
% OMI grid is bottom to top
      prs_lev=zeros(npxl,nscan,level);
      prs_lay=zeros(npxl,nscan,layer);
      for ipxl=1:npxl
         for iscan=1:nscan
            if(isnan(tm5_prs_sfc(ipxl,iscan,1)))
               continue
            end
            for k=1:layer
               prs_lev(ipxl,iscan,k)=tm5_prs_a(1,k)+tm5_prs_b(1,k)* ...
               tm5_prs_sfc(ipxl,iscan,1);
               prs_lay(ipxl,iscan,k)=(tm5_prs_a(1,k)+tm5_prs_a(2,k))/2.+ ...
               (tm5_prs_b(1,k)+tm5_prs_b(2,k))/2.*tm5_prs_sfc(ipxl,iscan,1);
            end
            prs_lev(ipxl,iscan,level)=tm5_prs_a(2,layer)+tm5_prs_b(2,layer)* ...
            tm5_prs_sfc(ipxl,iscan,1);
      end
% zenang(npxl,nscan,ntime) (deg)
      field='/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle';
      zenang=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% no2_slnt_col(npxl,nscan,ntime) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/scd_no2';
      no2_slnt_col=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% no2_slnt_col_err(npxl,nscan,ntime) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/scd_no2_uncertainty';
      no2_slnt_col_err=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% o3_slnt_col(npxl,nscan,ntime) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/scd_o3';
      o3_slnt_col=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% o3_slnt_col_err(npxl,nscan,ntime) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/scd_o3_uncertainty';
      o3_slnt_col_err=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% no2_strat_col(npxl,nscan,ntim) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/stratospheric_no2_vertical_column';
      no2_strat_col=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% no2_strat_col_err(npxl,nscan,ntim) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/stratospheric_no2_vertical_column_uncertainty';
      no2_strat_col_err=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% no2_total_col(npxl,nscan,ntim) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/total_no2_vertical_column';
      no2_total_col=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% no2_total_col_err(npxl,nscan,ntim) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/total_no2_vertical_column_uncertainty';
      no2_total_col_err=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% no2_summed_col(npxl,nscan,ntim) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/summed_no2_total_vertical_column';
      no2_summed_col=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% no2_summed_col_err(npxl,nscan,ntim) (molec/cm^2)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/summed_no2_total_vertical_column_uncertainty';
      no2_summed_col_err=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% amf_strat(npxl,nscan,ntime) ( )
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/amf_strat';
      amf_strat=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
%
% Loop through OMI data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,1995));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,1995));
      ocnt=0;
      icnt=0;
      for itim=1:ntim
         for iscan=1:nscan
            if(isnan(delta_time(iscan,itim)))
               continue
            end
            time_cur=ref_secs+delta_time(iscan,itim);
            [year,month,day,hour,minute,second]=invert_time_ref(time_cur,1995);
            yyyy_omi=year;
            mn_omi=month;
            dy_omi=day;
            hh_omi=hour;
            mm_omi=minute;
            ss_omi=second;
%            [jult]=convert_time_ref(year,month,day,hour,minute,second,1995);
            omidate=single(convert_time_ref(year,month,day,hour,minute,second,1995));
%
% Check time
            if(omidate<windate_min | omidate>windate_max)
%               fprintf('APM: min %d, omi %d, max %d \n',windate_min,omidate,windate_max);
               continue
            end
   	    for ipxl=1:npxl
%
% QA/AC
               if(isnan(proc_err_flag(ipxl,iscan,itim)) | proc_err_flag(ipxl,iscan,itim)~=0)
                  continue		
               end
               if(isnan(zenang(ipxl,iscan,itim)) | zenang(ipxl,iscan,itim)>=80.0)
                  continue
               end
               if(isnan(snow_ice_flag(ipxl,iscan,itim)) | (snow_ice_flag(ipxl,iscan,itim)>=10 & ...
               snow_ice_flag(ipxl,iscan,itim)~=255))
                  continue		
               end
               if(all(prs_lev(ipxl,iscan,:)==0.))
                  continue
               end
               if(isnan(amf_trop(ipxl,iscan,itim)) | isnan(amf_geom(ipxl,iscan,itim)) | ...
               amf_geom(ipxl,iscan,itim)==0)
                  continue		
               end
               if(amf_trop(ipxl,iscan,itim)/amf_geom(ipxl,iscan,itim)<=.2)
                  continue		
               end
               if(isnan(cld_rad_frac_no2(ipxl,iscan,itim)) | cld_rad_frac_no2(ipxl,iscan,itim)>0.5)
                  continue		
               end
               if(isnan(tm5_prs_sfc(ipxl,iscan,itim)))
                  continue
	        end
               if(any(isnan(prs_lay(ipxl,iscan,1:layer))))
                  continue
               end
               if(isnan(no2_trop_col(ipxl,iscan,itim)) | no2_trop_col(ipxl,iscan,itim)<=0)
                  continue
               end
               if(isnan(no2_slnt_col(ipxl,iscan,itim)) | no2_slnt_col(ipxl,iscan,itim)<=0)
                  continue
               end
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
	       xmdl_mn=lon_mdl(1,1);
	       if(xmdl_mn<0.)
	          xmdl_mn=xmdl_mn+360.;
               end
	       xmdl_mx=lon_mdl(nx_mdl,ny_mdl);
	       if(xmdl_mx<0.)
	          xmdl_mx=xmdl_mx+360.;
               end
%
               [xi,xj]=w3fb13(y_obser,x_obser,lat_mdl(1,1), ...
	       xmdl_mn,delx,cen_lon,truelat1,truelat2);
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
%                  fprintf('x_mdl_min, x_obs, x_mdl_max: %6.2f %6.2f %6.2f \n',xmdl_mn, ...
%                  x_obser,xmdl_mx)
%                  fprintf('y_mdl_min, y_obs, y_mdl_max: %6.2f %6.2f %6.2f \n',lat_mdl(1,1), ...
%                  y_obser,lat_mdl(nx_mdl,ny_mdl))
%                  fprintf('i_min %d j_min %d \n',i_min,j_min)
                  continue
               end
               if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
%                  fprintf('NO REJECT: i_min %d j_min %d \n',i_min,j_min)
                  continue
               end
%
% Save data to ascii file
               icnt=icnt+1;
               fprintf(fid,'OMI_NO2_Obs: %d %d %d \n',icnt,i_min,j_min);
               fprintf(fid,'%d %d %d %d %d %d \n',yyyy_omi, ...
	       mn_omi,dy_omi,hh_omi,mm_omi,ss_omi);
	       fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,iscan),lon(ipxl,iscan));
               fprintf(fid,'%d %d \n',layer,level);
               fprintf(fid,'%d \n',trop_indx(ipxl,iscan,itim));
               fprintf(fid,'%14.8g \n',amf_trop(ipxl,iscan,itim));
               fprintf(fid,'%14.8g \n',amf_total(ipxl,iscan,itim));
               fprintf(fid,'%14.8g %14.8g \n',no2_trop_col(ipxl,iscan,itim), ...
	       no2_trop_col_err(ipxl,iscan,itim));
               fprintf(fid,'%14.8g %14.8g \n',no2_slnt_col(ipxl,iscan,itim), ...
	       no2_slnt_col_err(ipxl,iscan,itim));
               fprintf(fid,'%14.8g ',avgk(1:layer,ipxl,iscan,itim))*amf_total(ipxl,iscan,itim);
               fprintf(fid,'\n');
 	       fprintf(fid,'%14.8g ',prs_lev(ipxl,iscan,1:level));
               fprintf(fid,'\n');
               if(all(prs_lev(ipxl,iscan,:)==0.))
                  fprintf('OMI PRS ZERO \n ')
                  fprintf('%14.8g ',prs_lev(ipxl,iscan,1:level))
                  fprintf(fid,'\n');
	          return
               end
	    end
         end
      end
   end
   clear nscan npxl ncnr ntim nply layer level nv
   clear lat lon time delta_time
   clear no2_trop_col no2_trop_col_err
   clear no2_slnt_col no2_slnt_col_err
   clear avgk amf_trop amf_total trop_indx
   clear tm5_prs_a tm5_prs_b tm5_sfc_prs
end

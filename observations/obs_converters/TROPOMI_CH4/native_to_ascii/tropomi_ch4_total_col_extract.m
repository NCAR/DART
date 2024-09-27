function tropomi_ch4_total_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
% col_avgk_lay(nlay,npxl,nscan) (none) (12, 215, 4172)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/column_averaging_kernel';
      col_avgk_lay=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
      temp_size=size(col_avgk_lay);
      nlay=temp_size(1);
      nlev=nlay+1;
      npxl=temp_size(2);
      nscan=temp_size(3);
% dofs(npxl,nscan) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/degrees_of_freedom';
      dofs=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% dofs_ch4(npxl,nscan) (none)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/degrees_of_freedom_methane';
      dofs_ch4=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% solar_zenith_angle(npxl,nscan) (degrees)
      field='/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle';
      zenang=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% alt_lev(nlev,npxl,nscan) (none) (13, 215, 4172) (top to bottom) (m)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/altitude_levels';
      alt_lev=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% alt_sfc(npxl,nscan) (none) (215, 4172)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude';
      alt_sfc=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% ch4_prof_prior(nlay,npxl,nscan) (none) (12, 215, 4172)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/methane_profile_apriori';
      ch4_prof_prior=double(ncread(file_in,field));
      fac_molec2cm2=ncreadatt(file_in,field,'multiplication_factor_to_convert_to_molecules_percm2');
      units=ncreadatt(file_in,field,'units');
% prs_sfc(npxl,nscan) (Pa)
      field='/PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure';
      prs_sfc=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% time_delta(nscan)
      field='/PRODUCT/delta_time';
      temp=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
      time_delt=double(temp)/1.e3;
% pixel(npxl)
      field='/PRODUCT/ground_pixel';
      temp=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');
% lat(npxl,nscan)
      field='/PRODUCT/latitude';
      lat=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% layer
      field='/PRODUCT/layer';
      temp=ncread(file_in,field);
% level
      field='/PRODUCT/level';
      temp=ncread(file_in,field);
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
% ch4_vmr(npxl,nscan)
      field='/PRODUCT/methane_mixing_ratio';
      ch4_vmr=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% ch4_vmr_bc(npxl,nscan)
      field='/PRODUCT/methane_mixing_ratio_bias_corrected';
      ch4_vmr_bc=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% ch4_vmr_err(npxl,nscan)
      field='/PRODUCT/methane_mixing_ratio_precision';
      ch4_vmr_err=double(ncread(file_in,field));
      units=ncreadatt(file_in,field,'units');
% qa_value(npxl,nscan)
      field='/PRODUCT/qa_value';
      qa_value=ncread(file_in,field); 
      fac_scale=ncreadatt(file_in,field,'scale_factor');   
      units=ncreadatt(file_in,field,'units');  
% scanline
      field='/PRODUCT/scanline';
      temp=ncread(file_in,field);
% time
      field='/PRODUCT/time';
      temp=ncread(file_in,field);
      units=ncreadatt(file_in,field,'units');  
      time=max(temp);
% time_utc
      field='/PRODUCT/time_utc';
      temp=ncread(file_in,field);
%
% Loop through TROPOMI data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2010));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2010));
      icnt=0;
      for ilin=1:nscan
         tropomidate=single(convert_time_ref(file_str_yy,file_str_mm,file_str_dd, ...
         0,0,0,2010))+time_delt(ilin);
         [yyyy_tropomi,mn_tropomi,dd_tropomi,hh_tropomi,mm_tropomi,ss_tropomi]= ...
	 invert_time_ref(tropomidate,2010);
%
% Check time
         if(tropomidate<windate_min | tropomidate>windate_max)
            continue
         end
         for ipxl=1:npxl
%
% QA/AC
            if(any(isnan(alt_lev(:,ipxl,ilin))) | any(alt_lev(:,ipxl,ilin)<=0))
               continue
            end
	    alt_lev_adj(:)=alt_lev(:,ipxl,ilin)-alt_sfc(ipxl,ilin);
%
            if (any(isnan(col_avgk_lay(:,ipxl,ilin))))
               continue
            end
%
            if (any(isnan(ch4_prof_prior(:,ipxl,ilin))) | any(ch4_prof_prior(:,ipxl,ilin)<=0))
               continue 
            end
%
            if(isnan(ch4_vmr(ipxl,ilin)) | ch4_vmr(ipxl,ilin)<=0)
               continue
            end
%
            if(isnan(ch4_vmr_err(ipxl,ilin)) | ch4_vmr_err(ipxl,ilin)<=0)
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
            xmdl_mx=lon_mdl(nx_mdl,ny_mdl);
            if(xmdl_mx<0.)
               xmdl_mx=xmdl_mx+360.;
            end
%
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
%               fprintf('FAILED DOMAIN TEST \n')	      
               continue
            end
            if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
%               fprintf('FAILED DOMAIN TEST \n')	      
               continue
            end
%
% Save data to ascii file
	    icnt=icnt+1;
	    fprintf('APM: Save the TROPOMI data point %d \n',icnt)
            fprintf(fid,'TROPOMI_CH4_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_tropomi, ...
            mn_tropomi,dd_tropomi,hh_tropomi,mm_tropomi,ss_tropomi);
            fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,ilin),lon(ipxl,ilin));
            fprintf(fid,'%d %d \n',nlay,nlev);
            fprintf(fid,'%14.8g %14.8g \n',alt_sfc(ipxl,ilin), ...
            prs_sfc(ipxl,ilin));
            fprintf(fid,'%14.8g ',alt_lev_adj(1:nlev));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',col_avgk_lay(1:nlay,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',ch4_prof_prior(1:nlay,ipxl,ilin));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g %14.8g \n',ch4_vmr(ipxl,ilin), ...
            ch4_vmr_err(ipxl,ilin));
         end
      end
   end
end

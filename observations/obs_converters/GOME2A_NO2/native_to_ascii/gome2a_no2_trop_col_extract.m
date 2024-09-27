function gome2a_no2_trop_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
      fprintf('READ GOME2A DATA \n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read GOME2A data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% get dimensions
      wfid=netcdf.open(file_in,'NC_NOWRITE');
      gfid=netcdf.inqNcid(wfid,'PRODUCT');
      dimid=netcdf.inqDimID(gfid,'scanline'); % 15571
      [name,nscan]=netcdf.inqDim(gfid,dimid); 
      dimid=netcdf.inqDimID(gfid,'ground_pixel'); % 1
      [name,npxl]=netcdf.inqDim(gfid,dimid); 
      dimid=netcdf.inqDimID(gfid,'corner'); % 4
      [name,ncnr]=netcdf.inqDim(gfid,dimid);
      dimid=netcdf.inqDimID(gfid,'gome2_across_track'); % 32
      [name,ntrk]=netcdf.inqDim(gfid,dimid);
      dimid=netcdf.inqDimID(gfid,'time'); % 1
      [name,ntim]=netcdf.inqDim(gfid,dimid);
      dimid=netcdf.inqDimID(gfid,'polynomial_exponents'); %5
      [name,nply]=netcdf.inqDim(gfid,dimid);
      dimid=netcdf.inqDimID(gfid,'layer'); % 34
      [name,layer]=netcdf.inqDim(gfid,dimid);
      level=layer+1;
      dimid=netcdf.inqDimID(gfid,'nv'); %2
      [name,nv]=netcdf.inqDim(gfid,dimid);
%
% lat (npxl,nscan)
      field='/PRODUCT/latitude';
      lat=ncread(file_in,field);
%
% lon (npxl,nscan)
      field='/PRODUCT/longitude';
      lon=ncread(file_in,field);
%
% time(ntim) (seconds since 1995-01-01 00:00:00)
      field='/PRODUCT/time';
      time=ncread(file_in,field);
%
% delta_time(nscan,ntim) (milliseconds from reference time of measurement)
      field='/PRODUCT/delta_time';
      delta_time=ncread(file_in,field);
      delta_time=delta_time/1000.;
      units=ncreadatt(file_in,field,'units');
%
% no2_vert_col_trop (npxl,nscan,ntim) (mollecules cm-2)
      field='/PRODUCT/tropospheric_no2_vertical_column';
      no2_vert_col_trop=ncread(file_in,field);
%
% no2_vert_col_err_trop (npxl,nscan,ntim)
      field='/PRODUCT/tropospheric_no2_vertical_column_uncertainty';
      no2_vert_col_err_trop=ncread(file_in,field);
%
% averaging_kernel (layer,npxl,nscan,ntim)
      field='/PRODUCT/averaging_kernel';
      avgk_lay=ncread(file_in,field);
%
% amf_trop (npxl,nscan,ntim)
      field='/PRODUCT/amf_trop';
      amf_trop=ncread(file_in,field);
%
% amf_total (npxl,nscan,ntim)
      field='/PRODUCT/amf_total';
      amf_total=ncread(file_in,field);
%
% trop_index (npxl,nscan,ntim)
      field='/PRODUCT/tm5_tropopause_layer_index';
      trop_index=ncread(file_in,field);
%
% tm5_prs_a (nv,layer) (hybrid pressure coefficient for Pa; 1 - lower grid bdy; 2 - upper grid bdy)
      field='/PRODUCT/tm5_pressure_level_a';
      tm5_prs_a=ncread(file_in,field);
      tm5_prs_a=tm5_prs_a/100.;
%
% tm5_prs_b (nv,layer) (hybrid pressure coefficient for Pa - lower grid bdy; 2 - upper grid bdy)
      field='/PRODUCT/tm5_pressure_level_b';
      tm5_prs_b=ncread(file_in,field);
%
% tm5_prs_sfc (npxl,nscan) (hPa)
      field='/PRODUCT/tm5_surface_pressure';
      tm5_prs_sfc=ncread(file_in,field);
%
% GOME2A grid is bottom to top
      prs_lev=zeros(npxl,nscan,level);
      prs_lay=zeros(npxl,nscan,layer);
      for ipxl=1:npxl
         for ilin=1:nscan
            if(isnan(tm5_prs_sfc(ipxl,ilin)))
               prs_lev(ipxl,ilin,:)=-99999.; 
               prs_lay(ipxl,ilin,:)=-99999.; 
               continue
            end
            for ilv=1:layer
               prs_lev(ipxl,ilin,ilv)=tm5_prs_a(1,ilv)+tm5_prs_b(1,ilv)* ...
               tm5_prs_sfc(ipxl,ilin);
               prs_lay(ipxl,ilin,ilv)=(tm5_prs_a(1,ilv)+tm5_prs_b(1,ilv)* ...
               tm5_prs_sfc(ipxl,ilin) + tm5_prs_a(2,ilv)+tm5_prs_b(2,ilv)* ...
               tm5_prs_sfc(ipxl,ilin))/2.;
            end
            prs_lev(ipxl,ilin,layer+1)=tm5_prs_a(2,layer)+tm5_prs_b(2,layer)* ...
            tm5_prs_sfc(ipxl,ilin);
         end
      end
%
% zenang (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle';
      zenang=ncread(file_in,field);
%
% no2_slnt_col (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/scd_no2';
      no2_slnt_col=ncread(file_in,field);
%
% no2_slnt_col_err (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/scd_no2_uncertainty';
      no2_slnt_col_err=ncread(file_in,field);
%
% o3_slnt_col (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/scd_o3';
      o3_slnt_col=ncread(file_in,field);
%
% o3_slnt_col_err (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/scd_o3_uncertainty';
      o3_slnt_col_err=ncread(file_in,field);
%
% no2_vert_col_strat; (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/stratospheric_no2_vertical_column';
      no2_vert_col_strat=ncread(file_in,field);
%
% no2_vert_col_err_strat (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/stratospheric_no2_vertical_column_uncertainty';
      no2_vert_col_err_total=ncread(file_in,field);
%
% no2_vert_col_total; (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/total_no2_vertical_column';
      no2_vert_col_total=ncread(file_in,field);
%
% no2_vert_col_err_total (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/total_no2_vertical_column_uncertainty';
      no2_vert_col_err_total=ncread(file_in,field);
%
% no2_vert_col_summed; (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/summed_no2_total_vertical_column';
      no2_vert_col_summedl=ncread(file_in,field);
%
% no2_vert_col_err_summed (npxl,nscan,nstep)
      field='/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/summed_no2_total_vertical_column_uncertainty';
      no2_vert_col_err_summed=ncread(file_in,field);
%
% Loop through GOME2A data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,1995));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,1995));
      icnt=0;
      for itim=1:ntim
         for iscan=1:nscan
            if(isnan(delta_time(iscan,itim)))
               continue
            end
            time_cur=ref_secs+delta_time(iscan,itim);
            [year,month,day,hour,minute,second]=invert_time_ref(time_cur,1995);
	    yyyy_gome2a=year;
	    mn_gome2a=month;
	    dy_gome2a=day;
	    hh_gome2a=hour;
	    mm_gome2a=minute;
	    ss_gome2a=second;
            gome2adate=single(convert_time_ref(year,month,day,hour,minute,second,1995));
%
% Check time
%            fprintf('APM: min %d, gome2a %d, max %d \n',windate_min,gome2adate,windate_max);
            if(gome2adate<windate_min | gome2adate>windate_max)
               continue
            end
	    for ipxl=1:npxl
%
% QA/AC
               if(any(isnan(prs_lev(itim,iscan,:))) | any(prs_lev(itim,iscan,:)<0))
                  continue
               end
%
               if(any(isnan(avgk_lay(:,ipxl,iscan,itim)))) 
                  continue
               end
%
               if(isnan(trop_index(ipxl,iscan,itim)) | trop_index(ipxl,iscan,itim)<=0)
                  continue
               end
%
               if(isnan(no2_vert_col_trop(ipxl,iscan,itim)) | no2_vert_col_trop(ipxl,iscan,itim)<=0)
                  continue
               end
%
               if(isnan(no2_vert_col_err_trop(ipxl,iscan,itim)) | no2_vert_col_err_trop(ipxl,iscan,itim)<=0)
                  continue
               end
%
               if(isnan(no2_vert_col_total(ipxl,iscan,itim)) | no2_vert_col_total(ipxl,iscan,itim)<=0)
                  continue
               end
%
               if(isnan(no2_vert_col_err_total(ipxl,iscan,itim)) | no2_vert_col_err_total(ipxl,iscan,itim)<=0)
                  continue
               end
%
               if(isnan(no2_slnt_col(ipxl,iscan,itim)) | no2_slnt_col(ipxl,iscan,itim)<=0)
                  continue
               end
%
               if(isnan(no2_slnt_col_err(ipxl,iscan,itim)) | no2_slnt_col_err(ipxl,iscan,itim)<=0)
                  continue
               end
%
%               if(isnan(o3_slnt_col(ipxl,iscan,itim)) | o3_slnt_col(ipxl,iscan,itim)<=0)
%                  continue
%               end
%
%               if(isnan(o3_slnt_col_err(ipxl,iscan,itim)) | o3_slnt_col_err(ipxl,iscan,itim)<=0)
%                  continue
%               end
%
               if(zenang(ipxl,iscan,itim)>=80.0)
                  continue
               end
%
% Check domain
% Input grid needs to be in degrees
% X coordinate is [0 to 360]
%
	       x_obser=lon(ipxl,iscan,itim);
               y_obser=lat(ipxl,iscan,itim);
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
%                  fprintf('DOMAIN ISSUE: i_min %d j_min %d \n',i_min,j_min)
                  continue
               end
               if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
%                  fprintf('NO REJECT: i_min %d j_min %d \n',i_min,j_min)
                  continue
               end
%
% Save data to ascii file
	       icnt=icnt+1;
               fprintf(fid,'GOME2A_NO2_Obs: %d %d %d \n',icnt,i_min,j_min);
               fprintf(fid,'%d %d %d %d %d %d \n',yyyy_gome2a, ...
               mn_gome2a,dy_gome2a,hh_gome2a,mm_gome2a,ss_gome2a);
               fprintf(fid,'%14.8f %14.8f \n',lat(ipxl,iscan,itim),lon(ipxl,iscan,itim));
               fprintf(fid,'%d %d \n',layer,level);
               fprintf(fid,'%d \n',trop_index(ipxl,iscan,itim));
               fprintf(fid,'%14.8g ',prs_lev(itim,iscan,1:level));
               fprintf(fid,'\n');
               fprintf(fid,'%14.8g ',avgk_lay(1:layer,ipxl,iscan,itim));
               fprintf(fid,'\n');
               fprintf(fid,'%14.8g %14.8g \n',no2_vert_col_trop(ipxl,iscan,itim), ...
               no2_vert_col_err_trop(ipxl,iscan,itim));
               fprintf(fid,'%14.8g %14.8g \n',no2_vert_col_total(ipxl,iscan,itim), ...
               no2_vert_col_err_total(ipxl,iscan,itim));
               fprintf(fid,'%14.8g \n',amf_trop(ipxl,iscan,itim));
               fprintf(fid,'%14.8g \n',amf_total(ipxl,iscan,itim));
               fprintf(fid,'%14.8g %14.8g \n',no2_slnt_col(ipxl,iscan,itim), ...
               no2_slnt_col_err(ipxl,iscan,itim));
               fprintf(fid,'%14.8g %14.8g \n',o3_slnt_col(ipxl,iscan,itim), ...
               o3_slnt_col_err(ipxl,iscan,itim));
            end
         end   
      end
   end
end

function mls_o3_profile_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)  
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
% Get the averaging kernel data
   file_avgk_path='/nobackupp11/amizzi/INPUT_DATA/FRAPPE_REAL_TIME_DATA/mls_o3_hdf_data/Averaging_Kernels_v4';
   file_avgk_1d0N='/MLS_v4_1D_AVK_0N.nc4';
   file_avgk_1d35N='/MLS_v4_1D_AVK_35N.nc4';
   file_avgk_1d70N='/MLS_v4_1D_AVK_70N.nc4';
%   file_avgk_2d0N='/MLS_v4_2D_AVK_0N.nc4';
%   file_avgk_2d35N='/MLS_v4_2D_AVK_35N.nc4';
%   file_avgk_2d70N='/MLS_v4_2D_AVK_70N.nc4';
%
% 1d_0N Data
   file_in=strcat(file_avgk_path,file_avgk_1d0N);
   wfid=netcdf.open(file_in,'NC_NOWRITE');
   gfid=netcdf.inqNcid(wfid,'O3');
   dimid=netcdf.inqDimID(gfid,'RetrievalLevel');
   [name,ret_nlay]=netcdf.inqDim(gfid,dimid); % 55
   dimid=netcdf.inqDimID(gfid,'TruthLevel');
   [name,tru_nlay]=netcdf.inqDim(gfid,dimid); % 55
   dimid=netcdf.inqDimID(gfid,'TruthPhi');
   [name,avk_nlat]=netcdf.inqDim(gfid,dimid);
%
% avgk_lay_1d0N(ret_nlay,tru_nlay)   
   field='/O3/avkv';
   avgk_lay_1d0N=ncread(file_in,field);
%
% ret_lay_0N(ret_nlay)
   field='/O3/RetrievalLevel';
   ret_lay_0N=ncread(file_in,field);
%
% tru_lay_0N(tru_nlay)
   field='/O3/TruthLevel';
   tru_lay_0N=ncread(file_in,field);
%
% avk_lat_0N(avk_nlat)
%   field='/O3/TruthPhi';
%   avk_lat_0N=ncread(file_in,field)
%
% 1d_35N Data
   file_in=strcat(file_avgk_path,file_avgk_1d35N);
%
% avgk_lay_1d35N(ret_nlay,tru_nlay)   
   field='/O3/avkv';
   avgk_lay_1d35N=ncread(file_in,field);
%
% ret_lay_35N(ret_nlay)
   field='/O3/RetrievalLevel';
   ret_lay_35N=ncread(file_in,field);
%
% tru_lay_35N(tru_nlay)
   field='/O3/TruthLevel';
   tru_lay_35N=ncread(file_in,field);
%
% avk_lat_35N(avk_nlat)
%   field='/O3/TruthPhi';
%   avk_lat_35N=ncread(file_in,field)
%
% 1d_70N Data
   file_in=strcat(file_avgk_path,file_avgk_1d70N);
%
% avgk_lay_1d70N(ret_nlay,tru_nlay)   
   field='/O3/avkv';
   avgk_lay_1d70N=ncread(file_in,field);
%
% ret_lay_70N(ret_nlay)
   field='/O3/RetrievalLevel';
   ret_lay_70N=ncread(file_in,field);
%
% tru_lay_70N(tru_nlay)
   field='/O3/TruthLevel';
   tru_lay_70N=ncread(file_in,field);
%
% avk_lat_70N(avk_nlat)
%   field='/O3/TruthPhi';
%   avk_lat_70N=ncread(file_in,field)
%
% Process satellite data
   for ifile=1:nfile
      clear time_start time_end
      clear nstep ntrk crnr layer level
      clear lat lon zenang time_utc
      clear o3_lay vert_col_total vert_col_trop
      clear prior_lay prior_err_lay prs_lev
      clear trop_indx avgk_lay avgk_diag_lay noise_corr info_content
      file_in=char(file_list(ifile));
      if(isempty(file_in))
         continue
      end
      indx=strfind(file_in,file_pre)-1;
      if(isempty(indx))
         continue
      end
      field='/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES';
      day=h5readatt(file_in,field,'GranuleDay');
      month=h5readatt(file_in,field,'GranuleMonth');
      year=h5readatt(file_in,field,'GranuleYear');
      time_start=h5readatt(file_in,field,'StartUTC');
      time_end=h5readatt(file_in,field,'EndUTC');
      tai93_0UTC=h5readatt(file_in,field,'TAI93At0zOfGranule');
%
      file_str_yy=str2double(time_start(1:4));
      file_str_mm=str2double(time_start(6:7));
      file_str_dd=str2double(time_start(9:10));
      file_str_hh=str2double(time_start(12:13));
      file_str_mn=str2double(time_start(15:16));
      file_str_ss=round(str2double(time_start(18:23)));
      file_end_yy=str2double(time_end(1:4));
      file_end_mm=str2double(time_end(6:7));
      file_end_dd=str2double(time_end(9:10));
      file_end_hh=str2double(time_end(12:13));
      file_end_mn=str2double(time_end(15:16));
      file_end_ss=round(str2double(time_end(18:23)));      
      file_str_secs=file_str_hh*3600 + file_str_mn*60 + file_str_ss;
      file_end_secs=file_end_hh*3600 + file_end_mn*60 + file_end_ss;
      fprintf('%d %s \n',ifile,file_in);
%      fprintf('file str secs %d cycle vend secs %d \n',file_secs,day_secs_end);
      if(file_str_secs>day_secs_end | file_end_secs<day_secs_beg)
%         fprintf('is %d <= %d, then process obs \n',file_str_secs,day_secs_end)
%         fprintf('is %d >= %d, then process obs \n',file_end_secs,day_secs_beg)
         continue
      end
      fprintf('READ MLS DATA \n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read MLS data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Time (ntime)
      field='/HDFEOS/SWATHS/O3/Geolocation Fields/Time';
      time=h5read(file_in,field);
      dmy=size(time);
      ntime=dmy(1);
%
% O3 Profile (layer,ntime) (VMR)     
      field='/HDFEOS/SWATHS/O3/Data Fields/L2gpValue';
      o3_obs=h5read(file_in,field);
      o3_units=h5readatt(file_in,field,'Units');
      clear dmy
      dmy=size(o3_obs);
      layer=dmy(1);
      level=layer+1;
%
% O3 Profile Error (layer,ntime)     
      field='/HDFEOS/SWATHS/O3/Data Fields/L2gpPrecision';
      o3_obs_err=h5read(file_in,field);
%
% O3 Profile Convergence (ntime)      
      field='/HDFEOS/SWATHS/O3/Data Fields/Convergence';
      o3_obs_conv=h5read(file_in,field);
%
% O3 Profile Quality (ntime)
      field='/HDFEOS/SWATHS/O3/Data Fields/Quality';
      o3_obs_qual=h5read(file_in,field);
%
% O3 Profile Status (ntime)
      field='/HDFEOS/SWATHS/O3/Data Fields/Status';
      o3_obs_status=h5read(file_in,field);
%
% lat (ntime)
      field='/HDFEOS/SWATHS/O3/Geolocation Fields/Latitude';
      lat=h5read(file_in,field);
%
% lon (ntime)
      field='/HDFEOS/SWATHS/O3/Geolocation Fields/Longitude';
      lon=h5read(file_in,field);
%
% prs_lay (layer) (hPa) Pressure is from bottom to top
      field='/HDFEOS/SWATHS/O3/Geolocation Fields/Pressure';
      prs_lay=h5read(file_in,field);
      units=h5readatt(file_in,field,'Units');
%
% zenang (ntime)
      field='/HDFEOS/SWATHS/O3/Geolocation Fields/SolarZenithAngle';
      zenang=h5read(file_in,field);
%
% O3 Column (ntime)
      field='/HDFEOS/SWATHS/O3 column/Data Fields/L2gpValue';
      o3_col_obs=h5read(file_in,field);
%
% O3 Column Error (ntime) 
      field='/HDFEOS/SWATHS/O3 column/Data Fields/L2gpPrecision';
      o3_col_obs_err=h5read(file_in,field);
%
% O3 Column Convergence (ntime)
      field='/HDFEOS/SWATHS/O3 column/Data Fields/Convergence';
      o3_col_obs_conv=h5read(file_in,field);
%
% O3 Column Quality (ntime)
      field='/HDFEOS/SWATHS/O3 column/Data Fields/Quality';
      o3_col_obs_qual=h5read(file_in,field);
%
% O3 Column Status (ntime)
      field='/HDFEOS/SWATHS/O3 column/Data Fields/Status';
      o3_col_obs_status=h5read(file_in,field);
%
% O3 A Priori (layer,ntime)
      field='/HDFEOS/SWATHS/O3-APriori/Data Fields/L2gpValue';
      o3_prior=h5read(file_in,field);
%
% O3 Prior Error (ntime)  
      field='/HDFEOS/SWATHS/O3-APriori/Data Fields/L2gpPrecision';
      o3_prior_err=h5read(file_in,field);
%
% O3 Prior Convergence (ntime)
      field='/HDFEOS/SWATHS/O3-APriori/Data Fields/Convergence';
      o3_prior_conv=h5read(file_in,field);
%
% O3 Prior Quality (ntime)
      field='/HDFEOS/SWATHS/O3-APriori/Data Fields/Quality';
      o3_prior_qual=h5read(file_in,field);
%
% O3 Prior Status (ntime)
      field='/HDFEOS/SWATHS/O3-APriori/Data Fields/Status';
      o3_prior_status=h5read(file_in,field);
%
% Loop through MLS data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2010));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2010));
      icnt=0;
      for itim=1:ntime
         curr_hrs=fix((time(itim)-tai93_0UTC)/3600);
         hrs_mod=mod((time(itim)-tai93_0UTC),3600);
         curr_min=fix(hrs_mod/60);
         curr_sec=mod(hrs_mod,60);
%
         yyyy_mls=year;
         mn_mls=month;
         dy_mls=day;
         hh_mls=curr_hrs;
         mm_mls=curr_min;
         ss_mls=round(curr_sec);
         mlsdate=single(convert_time_ref(yyyy_mls,mn_mls, ...
         dy_mls,hh_mls,mm_mls,ss_mls,2010));
         if(hh_mls==24 | mm_mls==60 | ss_mls==60)	 
            yyyy_tp=yyyy_mls;
            mn_tp=mn_mls;
            dy_tp=dy_mls;
            hh_tp=hh_mls;
            mm_tp=mm_mls;
            ss_tp=ss_mls;
            [yyyy_mls,mn_mls,dy_mls,hh_mls,mm_mls,ss_mls]= ...
            incr_time(yyyy_tp,mn_tp,dy_tp,hh_tp,mm_tp,ss_tp);
         end
%
% Check time
         if(mlsdate<windate_min | mlsdate>windate_max)
%            fprintf('min %d, mls %d, max %d \n',windate_min,mlsdate,windate_max)
            continue
         end
%
% QA/QC
%
         if(zenang(itim)>=80.0)
%            fprintf('zenang %6.2f \n',zenang(itim))
            continue
         end
%         if(isnan(o3_col_obs(itim)))
%            continue
%         end
%         for ilay=1:layer	 
%            if(isnan(o3_obs(ilay,itim)))
%	       fprinf('o3_obs is a NaN %d \n',o3_obs(ilay,itim)) 
%               continue
%            end
%         end
%	 
%         for ilay=1:layer	 
%            if(o3_obs(ilay,itim)<0.0)
%               fprintf('o3_obs is negative \n')
%               continue
%            end
%         end
%
% Check domain
% Input grid needs to be in degrees
% X coordinate is [0 to 360]
         x_obser=lon(itim);
         y_obser=lat(itim);
         if(x_obser<0.)
            x_obser=360.+x_obser;
         end
%
         xmdl_sw=lon_mdl(1,1);
         if(xmdl_sw<0.)
            xmdl_sw=xmdl_sw+360.;
         end
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
%            fprintf('x_mdl_min, x_obs, x_mdl_max: %6.2f %6.2f %6.2f \n',xmdl_sw, ...
%            x_obser,xmdl_mx)
%            fprintf('y_mdl_min, y_obs, y_mdl_max: %6.2f %6.2f %6.2f \n',lat_mdl(1,1), ...
%            y_obser,lat_mdl(nx_mdl,ny_mdl))
%            fprintf('i_min %d j_min %d \n',i_min,j_min)
            continue
         end
         if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
            fprintf('NO REJECT: i_min %d j_min %d \n',i_min,j_min)
            continue
         end
%
% Get the averaging kernel
         x_obser=lon(itim);
         y_obser=lat(itim);
         if(y_obser>=0 & y_obser<=35.)
           swt=35.-y_obser;
	   nwt=y_obser;
	   avgk_obs=(swt*avgk_lay_1d0N + nwt*avgk_lay_1d35N)/35.;
	   ret_lay_obs=(swt*ret_lay_0N + nwt*ret_lay_35N)/35.;
	   tru_lay_obs=(swt*tru_lay_0N + nwt*tru_lay_35N)/35.;
	 elseif (y_obser>35 & y_obser<=70.)
           swt=70.-y_obser;
	   nwt=y_obser-35.;
	   avgk_obs=(swt*avgk_lay_1d35N + nwt*avgk_lay_1d70N)/35.;
	   ret_lay_obs=(swt*ret_lay_35N + nwt*ret_lay_70N)/35.;
	   tru_lay_obs=(swt*tru_lay_35N + nwt*tru_lay_70N)/35.;
	 else
           fprintf('APM: Observation latitude is outside MLS AvgK latitude range \n')
	   continue
         end
%
% Save data to ascii file
         icnt=icnt+1;
         fprintf(fid,'MLS_O3_Obs: %d %d %d \n',icnt,i_min,j_min);
         fprintf(fid,'%d %d %d %d %d %d \n',yyyy_mls, ...
         mn_mls,dy_mls,hh_mls,mm_mls,ss_mls);
         fprintf(fid,'%14.8f %14.8f \n',lat(itim),lon(itim));
         fprintf(fid,'%d %d \n',layer,level);
         fprintf(fid,'%d %d \n',ret_nlay,tru_nlay);
%
% MLS prs_lay, ret_lay, and tru_lay are all the same
	 fprintf(fid,'%14.8g ',prs_lay(1:layer));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',ret_lay_obs(1:ret_nlay));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',tru_lay_obs(1:tru_nlay));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',o3_obs(1:layer,itim));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',o3_obs_err(1:layer,itim));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',o3_prior(1:layer,itim));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g \n',o3_prior_err(itim));
         for k=1:ret_nlay
            fprintf(fid,'%14.8g ',avgk_obs(k,1:tru_nlay));
            fprintf(fid,'\n');
	 end
%         fprintf(fid,'%14.8g \n',o3_col_obs(itim));
%         fprintf(fid,'%14.8g \n',o3_col_obs_err(itim));
      end
   end  
end

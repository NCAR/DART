function tempo_o3_cpsr_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
   command=strcat('ls'," ",'-1'," ",filein,'*.nc');
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
      clear time_start time_end
      clear nstep ntrk crnr layer level
      clear lat lon zenang time_utc
      clear o3_lay vert_col_total vert_col_trop
      clear prior_lay prior_err_lay prs_lev
      clear trop_indx avgk_lay noise_corr info_content
%
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
      file_str_secs=file_str_hh*3600 + file_str_mn*60 + file_str_ss;
      file_end_secs=file_end_hh*3600 + file_end_mn*60 + file_end_ss;
      fprintf('%d %s \n',ifile,file_in);
      fprintf('file str %d cycle end %d \n',file_str_secs,day_secs_end);
      fprintf('file end %d cycle str %d \n',file_end_secs,day_secs_beg);
%       
      if(file_str_secs>day_secs_end | file_end_secs<day_secs_beg)
         continue
      end
      fprintf('READ TEMPO DATA \n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read TEMPO data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% get dimensions
      wfid=netcdf.open(file_in,'NC_NOWRITE');
      [name,nstep]=netcdf.inqDim(wfid,0);
      [name,ntrk]=netcdf.inqDim(wfid,1);
      [name,crnr]=netcdf.inqDim(wfid,2);
      [name,layer]=netcdf.inqDim(wfid,3);
      [name,level]=netcdf.inqDim(wfid,4);
      netcdf.close(wfid);
%
% lat (ntrk,nstep)
      field='/geolocation/latitude';
      lat=ncread(file_in,field);
%
% lon (ntrk,nstep)
      field='/geolocation/longitude';
      lon=ncread(file_in,field);
%
% zenang (ntrk,nstep)
      field='/geolocation/solar_zenith_angle';
      zenang=ncread(file_in,field);
%
% time_utc (nstep) (seconds since 2001-01-01T12:00:00Z)
      field='/geolocation/time';
      time_utc=ncread(file_in,field);
%
% qa_value (ntrk,nstep)
%      field='/product/main_data_quality_flag';
%      qa_value=ncread(file_in,field);
%
% o3_lay (layer,ntrk,nstep)
      field='/product/ozone_profile';
      o3_lay=ncread(file_in,field);
      o3_lay=o3_lay*du2molpm2;
%
% vert_col_total (ntrk,nstep)
      field='/product/total_ozone_column';
      vert_col_total=ncread(file_in,field);
      vert_col_total=vert_col_total*du2molpm2;
%
% vert_col_trop (ntrk,nstep)
      field='/product/troposphere_ozone_column';
      vert_col_trop=ncread(file_in,field);
      vert_col_trop=vert_col_trop*du2molpm2;
%
% prior_lay (layer,ntrk,nstep)
      field='/support_data/ozone_apriori_profile';
      prior_lay=ncread(file_in,field);
      prior_lay=prior_lay*du2molpm2;
%
% prior_err_lay (layer,ntrk,nstep)
      field='/support_data/ozone_apriori_profile_error';
      prior_err_lay=ncread(file_in,field);
      prior_err_lay=prior_err_lay*du2molpm2;
%
% prs_lev (level,ntrk,nstep) (hPa)
      field='/support_data/ozone_profile_pressure';
      prs_lev=ncread(file_in,field);
%
% trop_index (ntrk,nstep)
      field='/support_data/tropopause_index';
      trop_index=ncread(file_in,field);
%
% avgk_lay (layer,layer,ntrk,nstep)
      field='/support_data/ozone_averaging_kernel';
      avgk_lay=ncread(file_in,field);
%
% noise_corr (layer,layer,ntrk,nstep)
      field='/support_data/ozone_noise_correlation_matrix';
      noise_corr=ncread(file_in,field);
      noise_corr=noise_corr*du2molpm2;
%
% info_content (ntrk,nstep)
      field='/support_data/ozone_information_content';
      info_content=ncread(file_in,field);
      fprintf('FINISH TEMPO DATA READ \n')
%
% Loop through TEMPO data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2000));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2000));
      icnt=0;
      for istep=1:nstep
         [jult]=jult_adjust_ref(time_utc(istep),2000,1,1,12,0,0);
         [yyyy_tempo,mn_tempo,dd_tempo,hh_tempo,mm_tempo,ss_tempo]=invert_time_ref(jult,2000);
         tempodate=single(convert_time_ref(yyyy_tempo,mn_tempo, ...
         dd_tempo,hh_tempo,mm_tempo,ss_tempo,2000));
%
% Check time
%	 fprintf('APM: time check %d %d %d \n',windate_min,tempodate,windate_max)
         if(tempodate<windate_min | tempodate>windate_max)
            continue
         end
         for ixtrk=1:ntrk
%
% QA/QC
            if(any(isnan(prs_lev(:,ixtrk,istep))) | any(prs_lev(:,ixtrk,istep)<0))
               continue
            end
%
            if(any(isnan(avgk_lay(:,:,ixtrk,istep))))
               continue
            end
%
            if(any(isnan(prior_lay(:,ixtrk,istep))) | any(prior_lay(:,ixtrk,istep)<0))
               continue
            end
%
            if(any(isnan(noise_corr(:,:,ixtrk,istep))))
               continue
            end
%
            if(isnan(vert_col_total(ixtrk,istep)) | vert_col_total(ixtrk,istep)<=0)
               continue
            end
%	 
            if(isnan(vert_col_trop(ixtrk,istep)) | vert_col_trop(ixtrk,istep)<=0)
               continue
            end
%	 
%	    if(qa_value(ixtrk,istep)~=0 | zenang(ixtrk,istep)>=80.0)
	    if(zenang(ixtrk,istep)>=80.0)
               continue
	    end
%
% Check domain
% Input grid needs to be in degrees
% X coordinate is [0 to 360]
	    x_obser=lon(ixtrk,istep);
            y_obser=lat(ixtrk,istep);
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
            fprintf(fid,'TEMPO_O3_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_tempo, ...
            mn_tempo,dd_tempo,hh_tempo,mm_tempo,ss_tempo);
            fprintf(fid,'%14.8f %14.8f \n',lat(ixtrk,istep),lon(ixtrk,istep));
            fprintf(fid,'%d %d \n',layer,level);
            fprintf(fid,'%d \n',trop_index(ixtrk,istep));
            fprintf(fid,'%14.8g ',prs_lev(1:level,ixtrk,istep));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',o3_lay(1:layer,ixtrk,istep));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',prior_lay(1:layer,ixtrk,istep));
            fprintf(fid,'\n');
	    for k=1:layer
               fprintf(fid,'%14.8g ',avgk_lay(k,1:layer,ixtrk,istep));
               fprintf(fid,'\n');
	    end
	    for k=1:layer
               fprintf(fid,'%14.8g ',noise_corr(k,1:layer,ixtrk,istep));
               fprintf(fid,'\n');
	    end
            fprintf(fid,'%14.8g \n',vert_col_trop(ixtrk,istep));
            fprintf(fid,'%14.8g \n',vert_col_total(ixtrk,istep));
	 end
      end
   end
end

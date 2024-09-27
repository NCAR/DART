function sciam_no2_total_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
      [name,nstep]=netcdf.inqDim(wfid,0);
      [name,ntrk]=netcdf.inqDim(wfid,1);
      [name,crnr]=netcdf.inqDim(wfid,2);
      [name,layer]=netcdf.inqDim(wfid,3);
      level=layer+1;
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
      field='/product/main_data_quality_flag';
      qa_value=ncread(file_in,field);
%
% no2_vert_col_trop (ntrk,nstep)
      field='/product/vertical_column_troposphere';
      no2_vert_col_trop=ncread(file_in,field);
%
% amf_total (ntrk,nstep)
      field='/support_data/amf_total';
      amf_total=ncread(file_in,field);
%
% amf_total_err (ntrk,nstep)
      field='/support_data/amf_total_uncertainty';
      amf_total_err=ncread(file_in,field);
%
% amf_trop (ntrk,nstep)
      field='/support_data/amf_troposphere';
      amf_trop=ncread(file_in,field);
%
% no2_slnt_col (ntrk,nstep)
      field='/support_data/fitted_slant_column';
      no2_slnt_col=ncread(file_in,field);
%
% no2_slnt_col_err (ntrk,nstep)
      field='/support_data/fitted_slant_column_uncertainty';
      no2_slnt_col_err=ncread(file_in,field);
%
% scat_wts (layer,ntrk,nstep)
      field='/support_data/scattering_weights';
      scat_wts=ncread(file_in,field);
%
% prs_sfc (ntrk,nstep)
      field='/support_data/surface_pressure';
      prs_sfc=ncread(file_in,field);
      Eta_A=ncreadatt(file_in,field,'Eta_A');
      Eta_B=ncreadatt(file_in,field,'Eta_B');
      for i=1:ntrk;
         for j=1:nstep;
            for k=1:level;
               prs_lev(k,i,j)=Eta_A(k)+(Eta_B(k)*prs_sfc(i,j));
            end
            for k=1:layer;
               prs_lay(k,i,j)=(prs_lev(k,i,j)+prs_lev(k+1,i,j))/2.;
            end
	 end
      end
%
% prs_trop (ntrk,nstep)
      field='/support_data/tropopause_pressure';
      prs_trop=ncread(file_in,field);
%
% no2_vert_col_total (ntrk,nstep)
      field='/support_data/vertical_column_total';
      no2_vert_col_total=ncread(file_in,field);
%
% no2_vert_col_total_err (ntrk,nstep)
      field='/support_data/vertical_column_total_uncertainty';
      no2_vert_col_total_err=ncread(file_in,field);
%
% Loop through GOME2A data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2000));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2000));
      icnt=0;
      for istep=1:nstep
	  [yyyy_gome2a,mn_gome2a,dy_gome2a,hh_gome2a,mm_gome2a,ss_gome2a]=invert_time_ref(time_utc(istep),2000);
%         if(int32(hh_gome2a)>23 | int32(mm_gome2a)>59 | ...
%         int32(ss_gome2a)>59)
%            [yyyy_gome2a,mn_gome2a,dy_gome2a,hh_gome2a, ...
%            mm_gome2a,ss_gome2a]=incr_time(yyyy_gome2a, ...
%      	 mn_gome2a,dy_gome2a,hh_gome2a,mm_gome2a,ss_gome2a);
%         end
%         fprintf('obs date/time %d %d %d %d %d %d \n',yyyy_gome2a, ...
%         mn_gome2a,dy_gome2a,hh_gome2a,mm_gome2a,ss_gome2a)
         gome2adate=single(convert_time_ref(yyyy_gome2a,mn_gome2a, ...
         dy_gome2a,hh_gome2a,mm_gome2a,ss_gome2a,2000));
%         fprintf('windate_min %d \n',windate_min)
%         fprintf('gome2a_dat %d \n',gome2adate)
%         fprintf('windate_max %d \n',windate_max)
%
% Check time
         if(gome2adate<windate_min | gome2adate>windate_max)
            continue
         end

	 for ixtrk=1:ntrk
%
% QA/AC
%	    if(qa_value(ixtrk,istep)~=0 | zenang(ixtrk,istep)>=80.0)
	    if(zenang(ixtrk,istep)>=80.0)
               continue
	    end
            if(isnan(no2_vert_col_trop(ixtrk,istep)) | no2_vert_col_trop(ixtrk,istep)<=0)
               continue
            end
%
% Check domain
% Input grid needs to be in degrees
% X coordinate is [0 to 360]
%		 
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
            fprintf(fid,'GOME2A_NO2_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_gome2a, ...
            mn_gome2a,dy_gome2a,hh_gome2a,mm_gome2a,ss_gome2a);
            fprintf(fid,'%14.8f %14.8f \n',lat(ixtrk,istep),lon(ixtrk,istep));
            fprintf(fid,'%d %d \n',layer,level);
            fprintf(fid,'%14.8g ',prs_lev(1:level,ixtrk,istep));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g ',scat_wts(1:layer,ixtrk,istep));
            fprintf(fid,'\n');
            fprintf(fid,'%14.8g %14.8g \n',no2_vert_col_total(ixtrk,istep), ...
            no2_vert_col_total_err(ixtrk,istep));
            fprintf(fid,'%14.8g %14.8g \n',no2_slnt_col(ixtrk,istep), ...
            no2_slnt_col_err(ixtrk,istep));
            fprintf(fid,'%14.8g %14.8g \n',amf_total(ixtrk,istep), ...
            amf_total_err(ixtrk,istep));
            fprintf(fid,'%14.8g \n',amf_trop(ixtrk,istep));
            fprintf(fid,'%14.8g \n',no2_vert_col_trop(ixtrk,istep));
            fprintf(fid,'%d \n',prs_trop(ixtrk,istep));         
         end
      end
   end
end

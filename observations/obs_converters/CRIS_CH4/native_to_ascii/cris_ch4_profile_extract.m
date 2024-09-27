function cris_ch4_profile_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)  
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
   command=strcat('ls'," ",'-1'," ",filein,'*')
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read CRIS data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
      fprintf('READ CRIS DATA \n')
      wfid=netcdf.open(file_in,'NC_NOWRITE');
      dimid=netcdf.inqDimID(wfid,'target'); % 39964
      [name,ntarg]=netcdf.inqDim(wfid,dimid); 
      dimid=netcdf.inqDimID(wfid,'level'); % 14
      [name,nlay]=netcdf.inqDim(wfid,dimid);
      nlev=nlay+1;
      dimid=netcdf.inqDimID(wfid,'datetime_utc_dim'); % 5
      [name,ndim]=netcdf.inqDim(wfid,dimid);
%
% lon (ntarg)
      field='/longitude';
      lon=ncread(file_in,field);
      for itarg=1:ntarg
         if(lon(itarg)<0.)
            lon(itarg)=360.+lon(itarg);
         end
      end
%
% lat (ntarg)
      field='/latitude';
      lat=ncread(file_in,field);
%
% date_time (ndim,ntarg) (YYYYMMDDHHMNSS)
      field='/datetime_utc';
      date_time=ncread(file_in,field);
%
% prs_lay (nlay,ntarg) (hPa)
      field='/pressure';
      prs_lay=ncread(file_in,field);
%
% ch4_lay (nlay,ntarg) (VMR)
      field='/x';
      ch4_lay=ncread(file_in,field);
%
% ch4_lay_prior (nlay,ntarg) (VMR)
      field='/observation_ops/xa';
      ch4_lay_prior=ncread(file_in,field);
%
% avgk_lay (nlay,nlay,ntarg) (ln(VMR_retr)/ln(VMR_true))
      field='/observation_ops/averaging_kernel';
      avgk_lay=ncread(file_in,field);
%
% err_cov_obs (nlay,nlay,ntarg) (ln(VMR_retr))
      field='/observation_ops/observation_error';
      err_cov_obs=ncread(file_in,field);
%
% dofs (ntarg) ( )
      field='/observation_ops/signal_dof';
      dofs=ncread(file_in,field);
%
% Loop through CRIS data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,2010));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,2010));
      icnt=0;
      for itarg=1:ntarg
         cris_yyyy=date_time(1,itarg);
         cris_mm=date_time(2,itarg);
         cris_dd=date_time(3,itarg);
         cris_hh=date_time(4,itarg);
         cris_mn=date_time(5,itarg);
         cris_ss=date_time(6,itarg);
         cris_time=single(convert_time_ref(cris_yyyy,cris_mm,cris_dd,cris_hh,cris_mn,cris_ss,2010));
%
% Check time
         if(cris_time<windate_min | cris_time>windate_max)
            fprintf('APM: min_date %d, cris_time %d, max_date %d \n',windate_min,cris_time,windate_max)
            continue
         end
%
% QA/QC
%
         reject=1
         for ilay=1:floor(nlay/2)
	    if(isnan(prs_lay(ilay,itarg)))
               continue
	    else
              reject=0
            end
         end
         if(reject==1)
            fprintf('APM: prs_lay is all NaNs \n')
            continue
         end
%      
	 for ilay=1:nlay
            if(isnan(ch4_lay(ilay,itarg)))
               fprintf('APM: ch4_lay has NaNs \n')
               continue
            end
         end
%
% Check domain
% Input grid needs to be in degrees
% X coordinate is [0 to 360]
         x_obser=lon(itarg);
         y_obser=lat(itarg);
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
            fprintf('x_mdl_min, x_obs, x_mdl_max: %6.2f %6.2f %6.2f \n',xmdl_sw, ...
            x_obser,xmdl_mx)
            fprintf('y_mdl_min, y_obs, y_mdl_max: %6.2f %6.2f %6.2f \n',lat_mdl(1,1), ...
            y_obser,lat_mdl(nx_mdl,ny_mdl))
            fprintf('i_min %d j_min %d \n',i_min,j_min)
            continue
         end
         if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
            fprintf('NO REJECT: i_min %d j_min %d \n',i_min,j_min)
            continue
         end
%
% Save data to ascii file
         icnt=icnt+1;
         fprintf(fid,'CRIS_CH4_Obs: %d %d %d \n',icnt,i_min,j_min);
         fprintf(fid,'%d %d %d %d %d %d \n',cris_yyyy, ...
         cris_mm,cris_dd,cris_hh,cris_mn,cris_ss);
         fprintf(fid,'%14.8f %14.8f \n',lat(itarg),lon(itarg));
         fprintf(fid,'%d %d \n',nlay,nlev);
         fprintf(fid,'%14.8g ',prs_lay(1:nlay,itarg));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',ch4_lay(1:nlay,itarg));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',ch4_lay_prior(1:nlay,itarg));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',dofs(itarg));
         fprintf(fid,'\n');
         for k=1:nlay
            fprintf(fid,'%14.8g ',avgk_lay(k,1:nlay,itarg));
            fprintf(fid,'\n');
	 end
         for k=1:nlay
            fprintf(fid,'%14.8g ',err_cov_obs(k,1:nlay,itarg));
            fprintf(fid,'\n');
	 end
      end
   end  
end

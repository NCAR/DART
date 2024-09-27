function modis_aod_total_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
   windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,1993));
   windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,1993));
%
   command=strcat('rm'," ",'-rf'," ",fileout);
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
%
% Convert DU to moles/m^2
   du2molpm2=4.4615e-4;
%
% Convert DU to molecules/m^2
%
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
% Process satellite data
   icnt=0;
   for ifile=1:nfile
      file_in=char(file_list(ifile));
      indx=strfind(file_in,file_pre)-1;
      if(isempty(indx))
         continue
      end
%
% time (TAI Time Seconds since 19930101 00:00:00.0 0)
      field='Scan_Start_Time';
      time=hdfread(file_in,field);
      tmp=size(time);
      nscan=tmp(1);
      npixl=tmp(2);
      if (time(1,1)>windate_max | time(nscan,npixl)<windate_min)
         continue
      end
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
      fprintf('READ MODIS DATA \n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read MODIS data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lat(nscan,npixl)
      field='Latitude';
      lat=hdfread(file_in,field);
% lon(nscan,npixl)
      field='Longitude';
      lon=hdfread(file_in,field);
% qa_lnd(nscan,npixl)
      field='Quality_Assurance_Land';
      qa_lnd=hdfread(file_in,field);
% qa_ocn(nscan,npixl)
      field='Quality_Assurance_Ocean';
      qa_ocn=hdfread(file_in,field);
% tau(nscan_npixl)
      field='Optical_Depth_Land_And_Ocean';
      tau=double(hdfread(file_in,field));
%      scale=h5readatt(file_in,field,'scale_factor');
      scale=.001;
      tau=tau*scale;
% taudb(nscan,npixl)
      field='Deep_Blue_Aerosol_Optical_Depth_550_Land';
      taudb=double(hdfread(file_in,field));
%      scale=h5readatt(file_in,field,'scale_factor');
      scale=.001;
      taudb=taudb*scale;
% qa_taudb(nxcan,npixl) (1 - Good; 2 - Very Good)
      field='Deep_Blue_Aerosol_Optical_Depth_550_Land_QA_Flag';
      qa_taudb=hdfread(file_in,field);
% taudb_best
      field='Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate';
      taudb_best=double(hdfread(file_in,field));
%      scale=h5readatt(file_in,field,'scale_factor');
      scale=.001;
      taudb_best=taudb_best*scale;
% taudb_err(nscan,npixl)
      field='Deep_Blue_Aerosol_Optical_Depth_550_Land_Estimated_Uncertainty';
      taudb_err=double(hdfread(file_in,field));
%      scale=h5readatt(file_in,field,'scale_factor');
      scale=.001;
      taudb_err=taudb_err*scale;
% zenang(nxcan,npixl)
      field='Solar_Zenith';
      zenang=hdfread(file_in,field);
%
% Loop through MODIS data
      for iscan=1:nscan
         for ipixl=1:npixl
%
% Check time
            if(time(iscan,ipixl)<windate_min | time(iscan,ipixl)>windate_max)
               continue
            end
%
% Check QA	 
            if(qa_taudb(iscan,ipixl)==0 )
               continue
            end
            if(isnan(taudb(iscan,ipixl)) | taudb(iscan,ipixl)<0)
               continue
            end
            if(isnan(taudb_best(iscan,ipixl)) | taudb_best(iscan,ipixl)<0)
               continue
            end
            if(isnan(taudb_err(iscan,ipixl)) | taudb_err(iscan,ipixl)<0)
               continue
            end
%
% Check domain
%
% Check domain
% Input grid needs to be in degrees
% X coordinate is [0 to 360]
%
            x_obser=lon(iscan,ipixl);
            y_obser=lat(iscan,ipixl);
            if(x_obser<0.)
               x_obser=360.+x_obser;
            end
%
            xmdl_sw=lon_mdl(1,1);
            if(xmdl_sw<0.)
               xmdl_sw=xmdl_sw+360.;
            end
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
% Invert MODIS time
            [yyyy_mod,mn_mod,dy_mod,hh_mod,mm_mod,ss_mod]=invert_time_ref(time(iscan,ipixl),1993);
%
% Save data to ascii file
            icnt=icnt+1;
            fprintf(fid,'MODIS_AOD_Obs: %d %d %d \n',icnt,i_min,j_min);
            fprintf(fid,'%d %d %d %d %d %d \n',yyyy_mod, ...
	    mn_mod,dy_mod,hh_mod,mm_mod,ss_mod);
	    fprintf(fid,'%14.8f %14.8f \n',lat(iscan,ipixl),lon(iscan,ipixl));
%	    fprintf('APM: taudb,best, err %14.8f %14.8f %14.8f \n',taudb(iscan,ipixl), ...
%	    taudb_best(iscan,ipixl),taudb_err(iscan,ipixl))
	    fprintf(fid,'%14.8f %14.8f %14.8f \n',taudb(iscan,ipixl),taudb_best(iscan,ipixl), ...
            taudb_err(iscan,ipixl));
         end
      end
      clear time lat lon qa_lnd qa_ocn zenang
      clear tau taudb qa_taudb taudb_best taudb_err 
   end
end

function mopitt_v5_co_profile_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
   layer=10
   prs_lev=[900. 800. 700. 600. 500. 400. 300. 200. 100. 50.];
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
      indx=strfind(file_in,file_pre)-1;
      if(isempty(indx))
         continue
      end
      year=str2double(file_in(indx+8:indx+11));
      month=str2double(file_in(indx+12:indx+13));
      day=str2double(file_in(indx+14:indx+15));
      file_hh=00;
      file_mm=00;
      file_secs_str=file_hh*60.*60. + file_mm*60.;
      file_hh=23;
      file_mm=59;
      file_secs_end=file_hh*60.*60. + file_mm*60. + 59.;

      fprintf('%d %s \n',ifile,file_in);
      fprintf('%d %d %d \n',day_secs_beg,file_secs_str,day_secs_end);
      fprintf('%d %d %d \n',day_secs_beg,file_secs_end,day_secs_end);
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read MOPITT data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Seconds in Day
      field='/Seconds in Day';
      temp=hdfread(file_in,field);
      secs_day=temp{:};
      temq=size(secs_day);
      numobs=temq(2);
% Latitude
      field='/Latitude';
      temp=hdfread(file_in,field);
      lat=temp{:};
% Longitude
      field='/Longitude';
      temp=hdfread(file_in,field);
      lon=temp{:};
% Solar Zenith Angle
      field='/Solar Zenith Angle';
      temp=hdfread(file_in,field);
      zen_ang=temp{:};
% Surface Pressure
      field='/Surface Pressure';
      temp=hdfread(file_in,field);
      prs_sfc=temp{:};
% Retrieved CO Total Column
      field='/Retrieved CO Total Column';
      temrr=hdfread(file_in,field);
      col_amt(:)=temrr(:,1);
      col_amt_err(:)=temrr(:,2);
% DOFs
      field='/Degrees of Freedom for Signal';
      temr=hdfread(file_in,field);
      dofs=temr{:};
% CO Mixing Ratio Profile (numobs,layerm,2)
      field='/Retrieved CO Mixing Ratio Profile';
      temrrr=hdfread(file_in,field);
      temo=size(temrrr);
      layerm=temo(2);
      co_retr_prf(:,:)=temrrr(:,:,1);
      co_retr_err_prf(:,:)=temrrr(:,:,2);
% CO Mixing Ratio Surface
      field='/Retrieved CO Surface Mixing Ratio';
      temrr=hdfread(file_in,field);
      co_retr_sfc=temrr(:,1);
      co_retr_err_sfc=temrr(:,2);
% Averaging Kernel(numobs,layer,layer)
      field='/Retrieval Averaging Kernel Matrix';
      temrrr=hdfread(file_in,field);
      avgk_lay(:,:,:)=temrrr(:,:,:);
% Prior Mixing Ratio Profile(numobs,layerm)
      field='/A Priori CO Mixing Ratio Profile';
      temrr=hdfread(file_in,field);
      co_prior_prf(:,:)=temrr(:,:,1);
      co_prior_err_prf(:,:)=temrr(:,:,2);
% Prior Mixing Ratio Surface
      field='/A Priori CO Surface Mixing Ratio';
      temr=hdfread(file_in,field);
      co_prior_sfc(:)=temr(:,1);
      co_prior_err_sfc(:)=temr(:,2);
% Retrieval Error Covariance(numobs,layer,layer)
      field='/Retrieval Error Covariance Matrix';
      temrrr=hdfread(file_in,field);
      cov_r(:,:,:)=temrrr(:,:,:);
%
% Loop through MOPITT data
      windate_min=single(convert_time(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn));
      windate_max=single(convert_time(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx));
      icnt=0;
%
      yyyy_mop=double(year);
      mn_mop=double(month);
      dy_mop=double(day);
      for iobs=1:numobs
         hh_mop=double(idivide(int32(secs_day(iobs)),3600));
         mm_mop=double(idivide(mod(int32(secs_day(iobs)),3600),60));
         ss_mop=double(int32(secs_day(iobs))-int32(hh_mop*3600+mm_mop*60));
         if(int32(hh_mop)>23 | int32(mm_mop)>59 | int32(ss_mop)>59)
            [yyyy_mop,mn_mop,dy_mop,hh_mop,mm_mop,ss_mop]=incr_time(yyyy_mop, ...
      	    mn_mop,dy_mop,hh_mop,mm_mop,ss_mop);
         end
         mopdate=single(convert_time(yyyy_mop,mn_mop,dy_mop,hh_mop,mm_mop,ss_mop));
%
% Check time
         if(mopdate<windate_min | mopdate>windate_max)
            continue
         end
%
% QA/AC
         if(isnan(col_amt(iobs)) | col_amt(iobs)<=0)
            continue
         end
%         if(isnan(prior_col_amt(iobs)) | prior_col_amt(iobs)<=0.)
%            continue
%         end
%
% Check domain
% Input grid needs to be in degrees
% X coordinate is [0 to 360]
%		 
         x_obser=lon(iobs);
         y_obser=lat(iobs);
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
%            fprintf('FAILED DOMAIN TEST \n')
         continue
         end
         if(i_min<1 | i_min>nx_mdl | j_min<1 | j_min>ny_mdl)
%           fprintf('FAILED DOMAIN TEST \n')
	    continue
	 end
%         fprintf('PASSED DOMAIN TEST \n')
%
% Save data to ascii file
         icnt=icnt+1;
         fprintf(fid,'MOPITT_CO_Profile_Obs: %d %d %d \n',icnt,i_min,j_min);
         fprintf(fid,'%d %d %d %d %d %d \n',yyyy_mop, ...
	 mn_mop,dy_mop,hh_mop,mm_mop,ss_mop);
	 fprintf(fid,'%14.8f %14.8f \n',lat(iobs),lon(iobs));
         fprintf(fid,'%d %d \n',layer,layer+1);
 	 fprintf(fid,'%14.8g \n',dofs(iobs));
 	 fprintf(fid,'%14.8g \n',prs_sfc(iobs));
         fprintf(fid,'%14.8g ',prs_lev(1:layer));
         fprintf(fid,'\n');
	 for k=1:layer
	    fprintf(fid,'%14.8g ',avgk_lay(iobs,k,1:layer));
            fprintf(fid,'\n');
	 end
	 fprintf(fid,'%14.8g \n',co_retr_sfc(iobs));
	 fprintf(fid,'%14.8g ',co_retr_prf(iobs,1:layerm));
         fprintf(fid,'\n');
	 fprintf(fid,'%14.8g \n',co_retr_err_sfc(iobs));
	 fprintf(fid,'%14.8g ',co_retr_err_prf(iobs,1:layerm));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g \n',co_prior_sfc(iobs));
         fprintf(fid,'%14.8g ',co_prior_prf(iobs,1:layerm));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g \n',co_prior_err_sfc(iobs));
	 fprintf(fid,'%14.8g ',co_prior_err_prf(iobs,1:layerm));
         fprintf(fid,'\n');
	 for k=1:layer
	    fprintf(fid,'%14.6g ',cov_r(iobs,k,1:layer));
            fprintf(fid,'\n');
         end
         fprintf(fid,'%14.8g \n',col_amt(iobs));
         fprintf(fid,'%14.8g \n',col_amt_err(iobs));
      end
      clear temp temq tempr temrr temrrr temo
      clear numobs layerm
      clear secs_day lat lon zen_ang col_amt col_amt_err dofs
      clear co_retr_prf co_retr_err_prf co_retr_sfc co_retr_err_sfc
      clear co_prior_prf co_prior_err_prf co_prior_sfc co_prior_err_sfc 
      clear avgk_lay cov_r
   end
end

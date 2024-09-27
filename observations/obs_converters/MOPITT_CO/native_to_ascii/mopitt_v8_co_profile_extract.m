function mopitt_v8_co_profile_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,path_mdl,file_mdl,cnx_mdl,cny_mdl)
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
   prss=[900. 800. 700. 600. 500. 400. 300. 200. 100.];
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
   for ifile=1:nfile
      file_in=char(file_list(ifile));
      indx=strfind(file_in,file_pre)-1;
      if(isempty(indx))
         continue
      end
      file_hh=00;
      file_mm=00;
      file_str_secs=file_hh*60.*60. + file_mm*60.;
      file_hh=23;
      file_mm=59;
      file_end_secs=file_hh*60.*60. + file_mm*60. + 59.;

      fprintf('%d %s \n',ifile,file_in);
      fprintf('file str %d cycle end %d \n',file_str_secs,day_secs_end);
      fprintf('file_end %d cycle str %d \n',file_end_secs,day_secs_beg);
%       
      if(file_str_secs>day_secs_end | file_end_secs<day_secs_beg)
         continue
      end
      fprintf('READ MOPITT DATA \n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read MOPITT data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For ntwo: 1 - value; 2 - uncertainty
% date data
      field='/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/';
      month=h5readatt(file_in,field,'Month');
      day=h5readatt(file_in,field,'Day');
      year=h5readatt(file_in,field,'Year');
      str_time=h5readatt(file_in,field,'StartDateTime');
% secs_day(ntim) 
      field='/HDFEOS/SWATHS/MOP02/Geolocation Fields/SecondsinDay';
      secs_day=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% time(ntim) (TAI time)
      field='/HDFEOS/SWATHS/MOP02/Geolocation Fields/Time';
      time=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% prior_prof_lay(ntwo,nlay,ntim) (ppbv)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/APrioriCOMixingRatioProfile';
      prior_prof_lay=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');
      arsize=size(prior_prof_lay);
      ntim=arsize(3);
      nlay=arsize(2);
      nlev=nlay+1;
      ntwo=arsize(1);
% prior_sfc(ntwo,ntim) (ppbv)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/APrioriCOSurfaceMixingRatio';
      prior_sfc=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% prior_col_amt(ntim) (ppbv) 
      field='/HDFEOS/SWATHS/MOP02/Data Fields/APrioriCOTotalColumn';
      prior_col_amt=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% avgk_row_sum_lev(nlev,ntim) (uses log10 VMR)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/AveragingKernelRowSums';
      avgk_row_sum_lev=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% dofs(ntim) (dimless)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/DegreesofFreedomforSignal';
      dofs=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% dry_col_amt(ntim) (molec/cm2)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/DryAirColumn';
      dry_col_amt=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% cov_m(nlev,nlev,ntim) 
      field='/HDFEOS/SWATHS/MOP02/Data Fields/MeasurementErrorCovarianceMatrix';
      cov_m=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% prs_grd(nlay) (hPa)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/PressureGrid';
      prs_grd=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% avgk_lev(nlev,nlev,ntim) (uses log10 VMR; ordered nrow,ncol,ntim)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/RetrievalAveragingKernelMatrix';
      avgk_lev=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% cov_r(nlev,nlev,ntim) 
      field='/HDFEOS/SWATHS/MOP02/Data Fields/RetrievalErrorCovarianceMatrix';
      cov_r=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% retr_prof_lay(ntwo,nlay,ntim) (ppbv)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/RetrievedCOMixingRatioProfile';
      retr_prof_lay=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% retr_sfc(ntwo,ntim) (ppbv)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/RetrievedCOSurfaceMixingRatio';
      retr_sfc=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% retr_col_amt(ntwo,ntim) (molec/cm2)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/RetrievedCOTotalColumn';
      retr_col_amt=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% cov_s(nlev,nlev,ntim) 
      field='/HDFEOS/SWATHS/MOP02/Data Fields/SmoothingErrorCovarianceMatrix';
      cov_s=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% zen_ang(ntim) (deg)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/SolarZenithAngle';
      zen_ang=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% prs_sfc(ntim) (hPa)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/SurfacePressure';
      prs_sfc=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% avgk_dim_col_lev(nlev,ntim) (uses log10 VMR)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/TotalColumnAveragingKernel';
      avgk_dim_col_lev=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% avgk_ndim_col_lev(nlev,ntim) (uses partial columns)
      field='/HDFEOS/SWATHS/MOP02/Data Fields/TotalColumnAveragingKernelDimless';
      avgk_ndim_col_lev=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% lat(ntim) (degrees)
      field='/HDFEOS/SWATHS/MOP02/Geolocation Fields/Latitude';
      lat=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
% lon(ntim) (degrees)
      field='/HDFEOS/SWATHS/MOP02/Geolocation Fields/Longitude';
      lon=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');
% prs_lay(nlay) (hPa) (900 hPa to 100 hPa)
      field='/HDFEOS/SWATHS/MOP02/Geolocation Fields/Pressure';
      prs_lay=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
%      nlay=size(prs_lay);
% prs_lev(nlev) (hPa) (1000 hPa to 100 hPa)
      field='/HDFEOS/SWATHS/MOP02/Geolocation Fields/Pressure2';
      prs_lev=h5read(file_in,field);
      units=h5readatt(file_in,field,'units');  
%      nlev=size(prs_lev);
%
% Loop through MOPITT data
      windate_min=single(convert_time_ref(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn,1993));
      windate_max=single(convert_time_ref(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx,1993));
      icnt=0;
      yyyy_mop=double(year);
      mn_mop=double(month);
      dy_mop=double(day);
      for itim=1:ntim
         hh_mop=double(idivide(int32(secs_day(itim)),3600));
         mm_mop=double(idivide(mod(int32(secs_day(itim)),3600),60));
         ss_mop=double(int32(secs_day(itim))-int32(hh_mop*3600+mm_mop*60));
         if(int32(hh_mop)>23 | int32(mm_mop)>59 | int32(ss_mop)>59)
            [yyyy_mop,mn_mop,dy_mop,hh_mop,mm_mop,ss_mop]=incr_time(yyyy_mop, ...
      	    mn_mop,dy_mop,hh_mop,mm_mop,ss_mop);
         end
         mopdate=single(convert_time_ref(yyyy_mop,mn_mop,dy_mop,hh_mop,mm_mop,ss_mop,1993));
%
% Check time
%         fprintf('min %d date %d max %d \n',windate_min,omidate,windate_max)
         if(mopdate<windate_min | mopdate>windate_max)
            continue
         end
%
% QA/AC
         if(isnan(retr_col_amt(1,itim)) | retr_col_amt(1,itim)<=0)
            continue
         end
         if(isnan(prior_col_amt(itim)) | prior_col_amt(itim)<=0.)
            continue
         end
         if all(prior_prof_lay(1,:,itim)==0)
	   continue
	 end
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
%
% Get corresponding indexes	 
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
% Fix pressure profile for surface pressure
	 prs_lev_sft=zeros(nlev);
	 prs_lev_sft(1:nlev)=prs_lev(1:nlev);
         for ilv=2:nlay
            if(ilv==1 || ilv==2 && prs_sfc(itim)>prs_lev_sft(ilv))
               kstr=1;
               nlev_sft=nlev;
               prs_lev_sft(1)=prs_sfc(itim);     
               break
            elseif(ilv==nlay && prs_sfc(itim)<prs_lev_sft(ilv+1))
               kstr=0;
               nlev_sft=0;
               break
            elseif(prs_sfc(itim)<=prs_lev_sft(ilv) && prs_sfc(itim)>prs_lev_sft(ilv+1))
               kstr=ilv;
               nlev_sft=nlev-kstr+1;
               prs_lev_sft(kstr)=prs_sfc(itim);     
               break	  
            end
         end
         if(nlev_sft==0)
            fprintf('APM: skip this ob - prs_sfc too low \n')
            continue
         end
	 nlay_sft=nlay-kstr+1;
         for ilv=1:kstr-1
            prs_lev_sft(ilv)=-999;
	 end
         for ilv=kstr:nlev
            prs_lev_sft(ilv-kstr+1)=prs_lev_sft(ilv);
         end
%
% Get shifted layer pressures
	 for ilv=1:nlay_sft 
	   prs_lay_sft(ilv)=(prs_lev_sft(ilv)+prs_lev_sft(ilv+1))/2.;
         end
%
% Shift remaining data	 
         retr_prof_sft=zeros(ntwo,nlay_sft);
         prior_prof_sft=zeros(ntwo,nlay_sft);
         avgk_sft=zeros(nlay_sft,nlay_sft);
         cov_m_sft=zeros(nlay_sft,nlay_sft);
         cov_r_sft=zeros(nlay_sft,nlay_sft);
         cov_s_sft=zeros(nlay_sft,nlay_sft);
	 for ilv=kstr:nlay
	    if(ilv-kstr+1==1) 
               retr_prof_sft(:,ilv-kstr+1)=retr_sfc(:,itim);
               prior_prof_sft(:,ilv-kstr+1)=prior_sfc(:,itim);
	    else
               retr_prof_sft(:,ilv-kstr+1)=retr_prof_lay(:,ilv,itim);
               prior_prof_sft(:,ilv-kstr+1)=prior_prof_lay(:,ilv,itim);
	    end
	 end
	 for ilv=kstr:nlay
	    for ilw=kstr:nlay
              avgk_sft(ilv-kstr+1,ilw-kstr+1)=avgk_lev(ilv+1,ilw+1,itim);
              cov_m_sft(ilv-kstr+1,ilw-kstr+1)=cov_m(ilv+1,ilw+1,itim);
              cov_r_sft(ilv-kstr+1,ilw-kstr+1)=cov_r(ilv+1,ilw+1,itim);
              cov_s_sft(ilv-kstr+1,ilw-kstr+1)=cov_s(ilv+1,ilw+1,itim);
            end
         end
%
% Save data to ascii file
         icnt=icnt+1;
         fprintf(fid,'MOPITT_CO_Obs: %d %d %d \n',icnt,i_min,j_min);
         fprintf(fid,'%d %d %d %d %d %d \n',yyyy_mop, ...
	 mn_mop,dy_mop,hh_mop,mm_mop,ss_mop);
	 fprintf(fid,'%14.8f %14.8f \n',lat(itim),lon(itim));
         fprintf(fid,'%d %d \n',nlay_sft,nlev_sft);
 	 fprintf(fid,'%14.8g \n',dofs(itim));
 	 fprintf(fid,'%14.8g \n',prs_sfc(itim));
         fprintf(fid,'%14.8g ',prs_lay_sft(1:nlay_sft));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',prs_lev_sft(1:nlev_sft));
         fprintf(fid,'\n');
	 for k=1:nlay_sft
            fprintf(fid,'%14.8g ',avgk_sft(k,1:nlay_sft));
            fprintf(fid,'\n');
	 end
	 fprintf(fid,'%14.8g ',retr_prof_sft(1,1:nlay_sft));
         fprintf(fid,'\n');
	 fprintf(fid,'%14.8g ',retr_prof_sft(2,1:nlay_sft));
         fprintf(fid,'\n');
         fprintf(fid,'%14.8g ',prior_prof_sft(1,1:nlay_sft));
         fprintf(fid,'\n');
	 fprintf(fid,'%14.8g ',prior_prof_sft(2,1:nlay_sft));
         fprintf(fid,'\n');
	 for k=1:nlay_sft
	    fprintf(fid,'%14.8g ',cov_s_sft(k,1:nlay_sft));
            fprintf(fid,'\n');
         end
	 for k=1:nlay_sft
            fprintf(fid,'%14.6g ',cov_r_sft(k,1:nlay_sft));
            fprintf(fid,'\n');
         end
	 for k=1:nlay_sft
            fprintf(fid,'%14.6g ',cov_m_sft(k,1:nlay_sft));
            fprintf(fid,'\n');
         end
         fprintf(fid,'%14.8g \n',retr_col_amt(1,itim));
         fprintf(fid,'%14.8g \n',retr_col_amt(2,itim));
         fprintf(fid,'%14.8g \n',prior_col_amt(itim));
      end
      clear prior_lay cld_prs col_amt rad_cld_frac avgk_lev 
      clear lat lon secs_day zen_ang time 
   end
end

function iasi_co_total_col_extract (filein,fileout,file_pre,cwyr_mn,cwmn_mn,cwdy_mn,cwhh_mn,cwmm_mn,cwss_mn,cwyr_mx,cwmn_mx,cwdy_mx,cwhh_mx,cwmm_mx,cwss_mx,clon_min,clon_max,clat_min,clat_max)
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
   lon_min=str2double(clon_min);
   lon_max=str2double(clon_max);
   lat_min=str2double(clat_min);
   lat_max=str2double(clat_max);
%
   command=strcat('rm'," ",'-rf'," ",fileout);
   [status]=system(command);
   fidot=fopen(fileout,'w');
%
   command=strcat('ls'," ",'-1'," ",filein,'*');
   [status,file_list_a]=system(command);
   file_list_b=split(file_list_a);
   file_list=squeeze(file_list_b)   
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
   fprintf('obs window str %d %d %d %d %d %d \n',wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn)
   fprintf('obs window end %d %d %d %d %d %d \n',wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx)
   fprintf('domain bounds %d %d %d %d \n',lat_min,lat_max,lon_min,lon_max)
%
   windate_min=single(convert_time(wyr_mn,wmn_mn,wdy_mn,whh_mn,wmm_mn,wss_mn));
   windate_max=single(convert_time(wyr_mx,wmn_mx,wdy_mx,whh_mx,wmm_mx,wss_mx));
%
   icnt=0;
   for ifile=1:nfile
      file_in=char(file_list(ifile));
      indx=strfind(file_in,file_pre)-1;
      if(isempty(indx))
         continue
      end
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
%      if((file_secs_str<day_secs_beg | file_secs>day_secs_end)
%         continue
%      end
%      fprintf('%d %s \n',ifile,file_in)
%      fprintf('%d %d %d \n',day_secs_beg,file_secs,day_secs_end)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read IASI data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
      fid = fopen(file_in,'r');
      iasi_data=textscan(fid,'%f %f %d %d %f %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
      fclose(fid);
      layer=19;
      level=20;
      nobs=size(iasi_data{1});
% lat
      lat=iasi_data{1,1};
% lon
      lon=iasi_data{1,2};
% date (yymmdd)
      date=iasi_data{1,3};
% time (hhmmss with no leading zeros)
      time=iasi_data{1,4};
% zen_ang
      zen_ang=iasi_data{1,5};
% dofs
      dofs=iasi_data{1,18};
% col_amt (molec/cm2)
      col_amt=iasi_data{1,21};  
% col_amt_err (molec/cm2)
      frac=iasi_data{1,22};
      for iobs=1:nobs
         col_amt_err(iobs)=frac(iobs)*col_amt(iobs);
      end
% prior_col_amt (molec/cm2) uses -999 when lower layers are missing
      for ilay=1:layer
	 prior_col_amt(:,ilay)=iasi_data{1,22+ilay}; 
      end
% avgk_lay (dimless) uses -999 when lower layers are missing
      for ilay=1:layer
	 avgk_lay(:,ilay)=iasi_data{1,41+ilay}; 
      end
%
      for iobs=1:nobs
	 year=floor(date(iobs)/10000);
	 month=floor(mod(date(iobs),10000)/100);
	 day=mod(date(iobs),100);
	 hour=floor(time(iobs)/10000);
	 minute=floor(mod(time(iobs),10000)/100);
	 second=mod(time(iobs),100);
	 secs_day(iobs)=hour*3600 + minute*60 + second;
	 yyyy_mop=double(year);
         mn_mop=double(month);
         dy_mop=double(day);
         hh_mop=double(hour);
	 mm_mop=double(minute);
	 ss_mop=double(second);
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
         if(isnan(prior_col_amt(iobs)) | prior_col_amt(iobs)<=0.)
            continue
         end
%
% Check domain
	 if(lon(iobs)<0)
	   lon(iobs)=lon(iobs)+360.;
	 end
	 if(lat(iobs)<lat_min | lat(iobs)>lat_max | ...
	 lon(iobs)<lon_min | lon(iobs)>lon_max)
            continue
         end
%
% Save data to ascii file
         icnt=icnt+1;
         fprintf(fidot,'IASI_CO_Obs: %d \n',icnt);
         fprintf(fidot,'%d %d %d %d %d %d \n',yyyy_mop, ...
	 mn_mop,dy_mop,hh_mop,mm_mop,ss_mop);
	 fprintf(fidot,'%14.8f %14.8f \n',lat(iobs),lon(iobs));
         fprintf(fidot,'%d %d \n',layer,level);
 	 fprintf(fidot,'%14.8g \n',dofs(iobs));
 	 fprintf(fidot,'%14.8g \n',col_amt(iobs));
 	 fprintf(fidot,'%14.8g \n',col_amt_err(iobs));
         fprintf(fidot,'%14.8g ',prior_col_amt(iobs,1:layer));
         fprintf(fidot,'\n');
         fprintf(fidot,'%14.8g ',avgk_lay(iobs,1:layer));
         fprintf(fidot,'\n');
      end
      clear prior_col_amt col_amt col_amt_err avgk_lay 
      clear lat lon secs_day zen_ang date time 
   end
end

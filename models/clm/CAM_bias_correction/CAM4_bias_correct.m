clear all
close all

%  2/10/23

%  CAM4 Reanalysis Product
%  V4 SIFMIP2 US-NR1 met forcing


% NR1
 SITE_lat=40.03;
 SITE_lon=-105.55+360; % degrees East
 SITE_doma_area=0.00109327562271889; % m^2 cell area (unique to grid location)
 SITE_doma_mask=1   ; % all values =1 

% Diagnostics ?  1=ON, 2=OFF
Diagnostics=1;

% CAM4:  1.875x2.5 degrees resolution every 6 hours for everything, 1997-2010
% CAM4 location:/glade/collections/rda/data/ds199.1/
% /glade/collections/rda/data/ds199.1/CAM_DATM.cpl_NINST.ha2x1dx6h.RUNYEAR.nc

% Site level extract CAM4:
% One file:
% noleap calendar
% SOLAR: a2x6h_Faxa_swvdf, a2x6h_Faxa_swndf, a2x6h_Faxa_swvdr, a2x6h_Faxa_swndr, a2x6h_Faxa_lwdn
% PRECIP: a2x6h_Faxa_rainc, a2x6h_Faxa_rainl, a2x6h_Faxa_snowc, a2x6h_Faxa_snowl
% MET:    a2x6h_Sa_pbot, a2x6h_Sa_shum, a2x6h_Sa_tbot, a2x6h_Sa_v, a2x6h_Sa_u
% unchanged: a2x6h_Sa_pslv, a2x6h_Sa_dens, a2x6h_Sa_ptem, a2x6h_Sa_z

path_CAM_forcing = '/glade/work/bmraczka/CAM4_NR1/';
% SIFMIP2 Tower met forcing files V2 had issues with NR1, use V4 files
% V2 files--> ver.2022.11.02.nc
% V4 files--> ver.2022.11.29.nc
path_towermet_forcing = '/glade/work/bmraczka/SIFMIP2/tower_met_forcing/'; 

% 'yearstr' is list of years where both CAM4 reanalysis and tower met forcing
% is available and the years we will create the bias-correct/scaled reanalysis product
yearstr={'2001','2002','2003','2004','2005','2006','2007','2008','2009','2010'};
% 'yeartower' adds previous year to cover the UTC to MST shift 
yeartower={'2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010'};


enstr={'0001','0002','0003','0004','0005','0006','0007','0008','0009','0010',...
         '0011','0012','0013','0014','0015','0016','0017','0018','0019','0020',...
         '0021','0022','0023','0024','0025','0026','0027','0028','0029','0030',...
         '0031','0032','0033','0034','0035','0036','0037','0038','0039','0040',...
         '0041','0042','0043','0044','0045','0046','0047','0048','0049','0050',...
         '0051','0052','0053','0054','0055','0056','0057','0058','0059','0060',...
         '0061','0062','0063','0064','0065','0066','0067','0068','0069','0070',...
         '0071','0072','0073','0074','0075','0076','0077','0078','0079','0080'};



% PLUMBER2 protocol time units --> LST (MST)


% PLUMBER2 protocol tower met forcing units
% Tair-->Kelvin
% Qair--> specific humidity (kg kg-1)
% Wind--> (m/s)
% SWdown --> (W/m2)
% LWdown ---> (W/m2)
% Precip --> (mm/s) or (kg m-2 s-1)
% Psurf  --> (Pa)


% Best to upload all years of tower forcing, then parse by year
TBOT_master=[];
SH_master=[];
WIND_master=[];
FSDS_master=[];
FLDS_master=[];
PRECTmms_master=[];
PSRF_master=[];
YEAR_master=[];

for ii = 1:length(yeartower);
%SIFMIP2 files
% V2 had issues with NR1, use V4 files
% V2 files--> ver.2022.11.02.nc
% V4 files--> ver.2022.11.29.nc
 
TBOT_dummy=squeeze(ncread([path_towermet_forcing 'Forcing_Tower_US-NR1_' yeartower{ii} '_ver.2022.11.29.nc'],'Tair')); 
SH_dummy=squeeze(ncread([path_towermet_forcing 'Forcing_Tower_US-NR1_' yeartower{ii} '_ver.2022.11.29.nc'],'Qair'));
WIND_dummy=squeeze(ncread([path_towermet_forcing 'Forcing_Tower_US-NR1_' yeartower{ii} '_ver.2022.11.29.nc'],'Wind'));
FSDS_dummy=squeeze(ncread([path_towermet_forcing 'Forcing_Tower_US-NR1_' yeartower{ii} '_ver.2022.11.29.nc'],'SWdown'));
FLDS_dummy=squeeze(ncread([path_towermet_forcing 'Forcing_Tower_US-NR1_' yeartower{ii} '_ver.2022.11.29.nc'],'LWdown'));
PRECTmms_dummy=squeeze(ncread([path_towermet_forcing 'Forcing_Tower_US-NR1_' yeartower{ii} '_ver.2022.11.29.nc'],'Precip'));
PSRF_dummy=squeeze(ncread([path_towermet_forcing 'Forcing_Tower_US-NR1_' yeartower{ii} '_ver.2022.11.29.nc'],'Psurf'));
YEAR_dummy=squeeze(ncread([path_towermet_forcing 'Forcing_Tower_US-NR1_' yeartower{ii} '_ver.2022.11.29.nc'],'year'));

TBOT_master=[TBOT_master;TBOT_dummy];
SH_master=[SH_master;SH_dummy];
WIND_master=[WIND_master;WIND_dummy];
FSDS_master=[FSDS_master;FSDS_dummy];
FLDS_master=[FLDS_master;FLDS_dummy];
PRECTmms_master=[PRECTmms_master;PRECTmms_dummy];
PSRF_master=[PSRF_master;PSRF_dummy];
YEAR_master=[YEAR_master;YEAR_dummy];

clear TBOT_dummy SH_dummy WIND_dummy FSDS_dummy FLDS_dummy PRECTmms_dummy PSRF_dummy YEAR_dummy

end

% Now identify each year to perform the correction
%  for ii = 1:length(yearstr);
   for ii = 5:5;   % TEST for a customized year. 

% Tower Forcing (PLUMBER):  LST (MST):  UTC-7
% CAM4  Forcing: UTC
% Tower Forcing must be advanced 7 hours to synchronize with CAM4 UTC

% Total shift of Tower Forcing: +7 hours

% CAM4 is 'noleap' calendar whereas tower forcing includes leap days
% Code accounts for this


% Select only indices of YEAR ii
indices=find(YEAR_master==str2num(yearstr{ii}));

% Pushing met forcing 7 hours forward (go back 14 indices)
ta_30=TBOT_master(indices(1)-14:indices(end)-14)';                         % TBOT (Kelvin)
q_30=SH_master(indices(1)-14:indices(end)-14)';                            % SH (kg kg-1)
wind_30=WIND_master(indices(1)-14:indices(end)-14)';                       % Total Wind (m/s) 
sw_30=FSDS_master(indices(1)-14:indices(end)-14)';                         % FSDS (W/m2)
lw_30=FLDS_master(indices(1)-14:indices(end)-14)';                         % FLDS (W/m2)
ppt_30=PRECTmms_master(indices(1)-14:indices(end)-14)';                    % PRECTmms  (mm/s) or (kg m-2 s-1)
ps_30=PSRF_master(indices(1)-14:indices(end)-14)';                         % PSRF (Pa)

% CAM4 is 'noleap' calendar, thus trim down tower met to 365 days to match
% Remove Feb 29th julian day: 60

if (str2num(yearstr{ii})==2000 | str2num(yearstr{ii})==2004 | str2num(yearstr{ii})==2008)  
leap_ind=59*48;
ta_30(leap_ind+1:leap_ind+48)=[];
q_30(leap_ind+1:leap_ind+48)=[];
wind_30(leap_ind+1:leap_ind+48)=[];
sw_30(leap_ind+1:leap_ind+48)=[];
lw_30(leap_ind+1:leap_ind+48)=[];
ppt_30(leap_ind+1:leap_ind+48)=[];
ps_30(leap_ind+1:leap_ind+48)=[];
end


clear indices

%% Generate 1,3,6 hourly averages to accomodate CAM4

ta_1hr=mean(reshape(ta_30,2,numel(ta_30)/2),1);        
q_1hr=mean(reshape(q_30,2,numel(q_30)/2),1);           
wind_1hr=mean(reshape(wind_30,2,numel(wind_30)/2),1);  
sw_1hr=mean(reshape(sw_30,2,numel(sw_30)/2),1);        
lw_1hr=mean(reshape(lw_30,2,numel(lw_30)/2),1);        
ppt_1hr=mean(reshape(ppt_30,2,numel(ppt_30)/2),1);     
ps_1hr=mean(reshape(ps_30,2,numel(ps_30)/2),1);


ta_6hr=mean(reshape(ta_1hr,6,numel(ta_1hr)/6),1);
q_6hr=mean(reshape(q_1hr,6,numel(q_1hr)/6),1);
wind_6hr=mean(reshape(wind_1hr,6,numel(wind_1hr)/6),1);
sw_6hr=mean(reshape(sw_1hr,6,numel(sw_1hr)/6),1);
lw_6hr=mean(reshape(lw_1hr,6,numel(lw_1hr)/6),1);
ppt_6hr=mean(reshape(ppt_1hr,6,numel(ppt_1hr)/6),1);
ps_6hr=mean(reshape(ps_1hr,6,numel(ps_1hr)/6),1);

clear ta_30 q_30 wind_30 sw_30 lw_30 ppt_30 ps_30

% CAM_DATM.cpl_NINST.ha2x1dx6h.RUNYEAR.nc
% SOLAR: a2x6h_Faxa_swvdf, a2x6h_Faxa_swndf, a2x6h_Faxa_swvdr, a2x6h_Faxa_swndr, a2x6h_Faxa_lwdn
% PRECIP: a2x6h_Faxa_rainc, a2x6h_Faxa_rainl, a2x6h_Faxa_snowc, a2x6h_Faxa_snowl
% MET:    a2x6h_Sa_pbot, a2x6h_Sa_shum, a2x6h_Sa_tbot, a2x6h_Sa_v, a2x6h_Sa_u
% Unchanged: a2x6h_Sa_pslv, a2x6h_Sa_dens, a2x6h_Sa_ptem, a2x6h_Sa_z
    
   for jj = 1:80;  %% 80  ensemble loop

       %Downloading all ensemble members for all  reanalysis variables
    
       Faxa_swndr(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Faxa_swndr'); 
       Faxa_swvdr(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Faxa_swvdr'); 
       Faxa_swndf(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Faxa_swndf');
       Faxa_swvdf(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Faxa_swvdf'); 
       Faxa_rainc(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Faxa_rainc'); 
       Faxa_rainl(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Faxa_rainl'); 
       Faxa_snowc(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Faxa_snowc'); 
       Faxa_snowl(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Faxa_snowl'); 
       Sa_u(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Sa_u'); 
       Sa_v(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Sa_v'); 
       Sa_tbot(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Sa_tbot'); 
       Sa_shum(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Sa_shum'); 
       Sa_pbot(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Sa_pbot');       
       Faxa_lwdn(jj,:)  = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Faxa_lwdn'); 
       Sa_z(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Sa_z');       
       Sa_ptem(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Sa_ptem'); 
       Sa_dens(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Sa_dens'); 
       Sa_pslv(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'],'a2x6h_Sa_pslv'); 

       

   end %% ensemble loop
     
      %% Initialize the ensemble mean scale values to the default of '0' mean no change
       Faxa_swndr_scale =zeros(1,length(Faxa_swndr(1,:))); 
       Faxa_swvdr_scale =zeros(1,length(Faxa_swvdr(1,:))); 
       Faxa_swndf_scale =zeros(1,length(Faxa_swndf(1,:))); 
       Faxa_swvdf_scale =zeros(1,length(Faxa_swvdf(1,:))); 

       Faxa_rainc_scale =zeros(1,length(Faxa_rainc(1,:))); 
       Faxa_rainl_scale =zeros(1,length(Faxa_rainl(1,:))); 
       Faxa_snowc_scale =zeros(1,length(Faxa_snowc(1,:))); 
       Faxa_snowl_scale =zeros(1,length(Faxa_snowl(1,:)));
       
       Faxa_rain_scale =zeros(1,length(Faxa_rainl(1,:)));
       Faxa_snow_scale =zeros(1,length(Faxa_snowl(1,:))); 

       Sa_u_scale =zeros(1,length(Sa_u(1,:))); 
       Sa_v_scale =zeros(1,length(Sa_v(1,:))); 

       Sa_tbot_scale =zeros(1,length(Sa_tbot(1,:))); 
       Sa_shum_scale =zeros(1,length(Sa_shum(1,:))); 
       Sa_pbot_scale =zeros(1,length(Sa_pbot(1,:))); 
       Faxa_lwdn_scale  =zeros(1,length(Faxa_lwdn(1,:))); 

       Sa_z_scale =zeros(1,length(Sa_z(1,:))); 
       Sa_ptem_scale =zeros(1,length(Sa_ptem(1,:))); 
       Sa_dens_scale =zeros(1,length(Sa_dens(1,:))); 
       Sa_pslv_scale =zeros(1,length(Sa_pslv(1,:)));
      
       % Calculate the ensemble mean for each CAM4 variable
        Faxa_swndr_mean= mean(Faxa_swndr,1);
        Faxa_swvdr_mean= mean(Faxa_swvdr,1);
        Faxa_swndf_mean= mean(Faxa_swndf,1);
        Faxa_swvdf_mean= mean(Faxa_swvdf,1);
     
        Faxa_rainc_mean= mean(Faxa_rainc,1);
        Faxa_rainl_mean= mean(Faxa_rainl,1);
        Faxa_snowc_mean= mean(Faxa_snowc,1);
        Faxa_snowl_mean= mean(Faxa_snowl,1);

        Sa_u_mean=mean(Sa_u,1);
        Sa_v_mean=mean(Sa_v,1);

        % Calculate CAM4 total wind vector for comparison with tower forcing
        Sa_wind_mean= (Sa_v_mean.^2+Sa_u_mean.^2).^0.5;

       Sa_tbot_mean = mean(Sa_tbot,1);
       Sa_shum_mean = mean(Sa_shum,1);
       Sa_pbot_mean = mean(Sa_pbot,1);
       Faxa_lwdn_mean = mean(Faxa_lwdn,1);


     % Need this initialization for solar loop
       Faxa_rain=Faxa_rainl;
       Faxa_snow=Faxa_snowl;


      %% Calculate the scaled corrections for each year 

      %% SOLAR 6 hour resolution

      % Identify year of interest for tower met forcing      

      for ind=1:(length(sw_6hr));

        if (Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind)+Faxa_swndf_mean(ind)+Faxa_swvdr_mean(ind)>0)      
      
           Faxa_swndr_scale(ind)= Faxa_swndr_mean(ind) .* (sw_6hr(ind)./(Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind) ...
                                                             +Faxa_swndf_mean(ind)+Faxa_swvdf_mean(ind)));
           Faxa_swvdr_scale(ind)= Faxa_swvdr_mean(ind) .* (sw_6hr(ind)./(Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind) ...
                                                             +Faxa_swndf_mean(ind)+Faxa_swvdf_mean(ind)));
           Faxa_swndf_scale(ind)= Faxa_swndf_mean(ind) .* (sw_6hr(ind)./(Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind) ...
                                                             +Faxa_swndf_mean(ind)+Faxa_swvdf_mean(ind)));
           Faxa_swvdf_scale(ind)= Faxa_swvdf_mean(ind) .* (sw_6hr(ind)./(Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind) ...
                                                             +Faxa_swndf_mean(ind)+Faxa_swvdf_mean(ind)));
        else
        % Do nothing, default scale = 0
        end
      end 
      %% PRECIP 6  hour resolution
      
      
      for ind=1:(length(ppt_6hr));
      
      % Precip loop needs to query temperature to re-assign to snow or rain

          if ta_6hr(ind)>273.15
             % Forcing all unadjusted precip ensemble members to new variable Faxa_rain according to adjusted temperature
             Faxa_rain(:,ind)=Faxa_rainc(:,ind)+Faxa_rainl(:,ind)+Faxa_snowc(:,ind)+Faxa_snowl(:,ind);
             Faxa_snow(:,ind)=0;

             % Adjustment will only be applied to total rain
             Faxa_rainc_scale(ind)= 0;
             Faxa_rainl_scale(ind)= 0;
             Faxa_snowc_scale(ind)= 0;
             Faxa_snowl_scale(ind)= 0;
             Faxa_rain_scale(ind)= ppt_6hr(ind);
          else
            % Forcing all unadjusted precip ensemble members to new variable Faxa_ snow according to adjusted temperature
            Faxa_snow(:,ind)= Faxa_rainc(:,ind)+Faxa_rainl(:,ind)+Faxa_snowc(:,ind)+Faxa_snowl(:,ind);
            Faxa_rain(:,ind)=0;

            % Adjustment will only be applied to total snow
             Faxa_rainc_scale(ind)=0;
             Faxa_rainl_scale(ind)=0;
             Faxa_snowc_scale(ind)=0;
             Faxa_snowl_scale(ind)=0;
             Faxa_snow_scale(ind)= ppt_6hr(ind);
          end

      end

       Faxa_rain_mean= mean(Faxa_rain,1);
       Faxa_snow_mean= mean(Faxa_snow,1);

       % Remaining 6 hour variables

       Sa_u_scale = Sa_u_mean.*(wind_6hr./Sa_wind_mean); 
       Sa_v_scale = Sa_v_mean.*(wind_6hr./Sa_wind_mean);
       Sa_tbot_scale = Sa_tbot_mean.*(ta_6hr./Sa_tbot_mean);
       Sa_shum_scale = Sa_shum_mean.*(q_6hr./Sa_shum_mean);
       Sa_pbot_scale = Sa_pbot_mean.*(ps_6hr./Sa_pbot_mean);
       Faxa_lwdn_scale = Faxa_lwdn_mean.*(lw_6hr./Faxa_lwdn_mean);

       %(Unchanged variables - no scaling)
       %Sa_z_scale    =  Sa_z_mean;
       %Sa_ptem_scale =  Sa_ptem_mean;
       %Sa_dens_scale =  Sa_desn_mean;
       %Sa_pslv_scale =  Sa_pslv_mean;
       %Sa_topo_scale =  Sa_topo_mean;

   %% Apply the scaled corrections identically to each ensemble member, then generate updated netcdf
   %% ensemble_scaled=ensemble + (ensemble_mean_scaled-ensemble_mean)

   %% Put clamps on non-physical values.

       % SOLAR
       % e.g Faxa_swndr(ensemble(80),time(6hr))
       Faxa_swndr_adjust  = Faxa_swndr+repmat(Faxa_swndr_scale-Faxa_swndr_mean,80,1); Faxa_swndr_adjust(Faxa_swndr_adjust<0)=0;  
       Faxa_swvdr_adjust  = Faxa_swvdr+repmat(Faxa_swvdr_scale-Faxa_swvdr_mean,80,1); Faxa_swvdr_adjust(Faxa_swvdr_adjust<0)=0;
       Faxa_swndf_adjust  = Faxa_swndf+repmat(Faxa_swndf_scale-Faxa_swndf_mean,80,1); Faxa_swndf_adjust(Faxa_swndf_adjust<0)=0;
       Faxa_swvdf_adjust  = Faxa_swvdf+repmat(Faxa_swvdf_scale-Faxa_swvdf_mean,80,1); Faxa_swvdf_adjust(Faxa_swvdf_adjust<0)=0;

       %PRECIP
       % Purposely setting convective rain and snow to zero
       % Purposely setting snow and rain adjusted variables to large scale snow (snowl) and rain (rainl)
       % These will carry the ensemble spread for snow/rain
 
       Faxa_rainc_adjust = Faxa_rainc+repmat(Faxa_rainc_scale-Faxa_rainc_mean,80,1); Faxa_rainc_adjust(:,:)=0;
       Faxa_rainl_adjust = Faxa_rain+repmat(Faxa_rain_scale-Faxa_rain_mean,80,1); Faxa_rainl_adjust(Faxa_rainl_adjust<0)=0;
       Faxa_snowc_adjust = Faxa_snowc+repmat(Faxa_snowc_scale-Faxa_snowc_mean,80,1); Faxa_snowc_adjust(:,:)=0;
       Faxa_snowl_adjust = Faxa_snow+repmat(Faxa_snow_scale-Faxa_snow_mean,80,1); Faxa_snowl_adjust(Faxa_snowl_adjust<0)=0;

       % meridional and zonal winds allowed to go negative
       Sa_u_adjust = Sa_u+repmat(Sa_u_scale-Sa_u_mean,80,1);
       Sa_v_adjust = Sa_v+repmat(Sa_v_scale-Sa_v_mean,80,1); 

       % Remaining variables
       Sa_tbot_adjust = Sa_tbot+repmat(Sa_tbot_scale-Sa_tbot_mean,80,1); Sa_tbot_adjust(Sa_tbot_adjust<0)=0;
       Sa_shum_adjust = Sa_shum+repmat(Sa_shum_scale-Sa_shum_mean,80,1); Sa_shum_adjust(Sa_shum_adjust<0)=0;
       Sa_pbot_adjust = Sa_pbot+repmat(Sa_pbot_scale-Sa_pbot_mean,80,1); Sa_pbot_adjust(Sa_pbot_adjust<0)=0;
       Faxa_lwdn_adjust = Faxa_lwdn+repmat(Faxa_lwdn_scale-Faxa_lwdn_mean,80,1); Faxa_lwdn_adjust(Faxa_lwdn_adjust<0)=0;

       % (3hr STATE Unchanged)
       %Sa_z = Sa_z+repmat(Sa_z_scale-Sa_z_mean,80,1); 
       %Sa_ptem = Sa_ptem+repmat(Sa_ptem_scale-Sa_ptem_mean,80,1);
       %Sa_dens = Sa_dens+repmat(Sa_dens_scale-Sa_dens_mean,80,1);
       %Sa_pslv = Sa_pslv+repmat(Sa_pslv_scale-Sa_pslv_mean,80,1);
       %Sa_topo = Sa_topo+repmat(Sa_topo_scale-Sa_topo_mean,80,1);

if Diagnostics==1

  % Diagnostic Figures Only
     % Comparing original CAM4 and Tower met forcing
     % Check for timing synchronization, variable magnitude checks 
     
       figure(1)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       j=plot(repelem(mean(Faxa_swndr,1),1,6)'+ ...
              repelem(mean(Faxa_swndf,1),1,6)'+ ...
              repelem(mean(Faxa_swvdr,1),1,6)'+ ...
              repelem(mean(Faxa_swvdf,1),1,6)', ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(sw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Total Shortwave (W m^{-2})','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr{ii}] , 'FontSize',fontpt);

       axes(ha(2));
       j=plot(repelem(mean(Faxa_rainl,1),1,6)'+ ...
              repelem(mean(Faxa_rainc,1),1,6)'+ ...
              repelem(mean(Faxa_snowl,1),1,6)'+ ...
              repelem(mean(Faxa_snowc,1),1,6)', ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ppt_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)

       axes(ha(3));
       j=plot(repelem(mean(Sa_tbot,1),1,6), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ta_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Temperature (K)','FontSize', fontpt)

       axes(ha(4));
       j=plot(repelem(mean(Sa_shum,1),1,6), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(q_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Specific Humidity (kg kg^{-1})','FontSize', fontpt) 
       linkaxes(ha(1:4), 'x'); 
      
       figure(2)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(3,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       j=plot(repelem(mean(Sa_u,1),1,6)', ...
              '--', 'color',rgb('darksalmon'), 'linewidth', mean_width); hold on
       k=plot(repelem(mean(Sa_v,1),1,6)', ...
              '-', 'color',rgb('darksalmon'), 'linewidth', mean_width);
       l=plot((repelem(mean(Sa_u,1),1,6)'.^2+ ...
               repelem(mean(Sa_v,1),1,6)'.^2).^0.5, ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width);

       m=plot(wind_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k,l,m], 'Zonal, CAM4 Ens Mean', 'Meridional, CAM4 Ens Mean', ...
                     'Total, CAM4 Ens Mean','Tower Met','FontSize',fontpt);

       ylabel('Wind (m s^{-1})','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr{ii}] , 'FontSize',fontpt);


       axes(ha(2));
       j=plot(repelem(mean(Sa_tbot,1),1,6), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ta_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Temperature (K)','FontSize', fontpt)

       axes(ha(3));
       j=plot(repelem(mean(Faxa_rainl,1),1,6)'+ ...
              repelem(mean(Faxa_rainc,1),1,6)'+ ...
              repelem(mean(Faxa_snowl,1),1,6)'+ ...
              repelem(mean(Faxa_snowc,1),1,6)', ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ppt_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)

       linkaxes(ha(1:3), 'x');


       figure(3)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(2,1, [.02,.08],[.10,.06],[.08,.06]);

       axes(ha(1));
       j=plot(repelem(mean(Sa_pbot,1),1,6), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ps_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Pressure (Pa)','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr{ii}] , 'FontSize',fontpt);

       axes(ha(2));
       j=plot(repelem(mean(Faxa_lwdn,1),1,6), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(lw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('LW radiation (W m^{-1})','FontSize', fontpt)
       linkaxes(ha(1:2), 'x');

display('  ')
display('Finished tower CAM4 diagnostic plots. Press enter to proceed to adjusted CAM4 plots..')
display('  ')
pause
close all

     % Compare CAM4 original vs CAM4 scaled, sanity check that algorithm is behaving as expected

     % 1)swndr 2)swvdr 3)swndf 4)swvdf 5)rainc 6)rainl 7)snowc 8)snowl  9)u 10)v
     % 11)tbot    12)shum      13)pbot    14)lwdn      
     % CAM4 mean (red); CAM4 scaled mean (green);    
     % CAM4 ensemble(light red) ; CAM4 scaled ensemble (light green)
     % Figure 1 Summer Month
     % Figure 2 Winter Month

       
     % Matrices (80, 2920) --> (80, 8760); 
     % repelem(Faxa_rainc_adjust,1,3);       

       figure(4)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3)/2 scrsz(4)]);
       fontpt=16;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       ax1=plot(Faxa_swndr(:,188*4:195*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndr_adjust(:,188*4:195*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width); 
       j=plot(mean(Faxa_swndr(:,188*4:195*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); 
       k=plot(mean(Faxa_swndr_adjust(:,188*4:195*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on 
       legend( [j,k], 'CAM4 Mean', 'Adjusted CAM4 Mean','FontSize',fontpt);       
       ylabel('Near-IR, Direct (W  m^{-2})','FontSize', fontpt)
       title(['CAM4 scaling against tower NR1 Met: July ', yearstr{ii}] , 'FontSize',fontpt);
       %xticks([24:24:24*18]);
       set(ha(1),'xticklabel',[]) 

       axes(ha(2));

       ax2=plot(Faxa_rain(:,188*4:195*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_rainl_adjust(:,188*4:195*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width); 
       plot(mean(Faxa_rain(:,188*4:195*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_rainl_adjust(:,188*4:195*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Rain (mm  s^{-1})','FontSize', fontpt)
       %xticks([24:24:24*18]);
       set(ha(2),'xticklabel',[])

       axes(ha(3));

       ax3=plot(Sa_tbot(:,188*4:195*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_tbot_adjust(:,188*4:195*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width); 
       plot(mean(Sa_tbot(:,188*4:195*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); 
       plot(mean(Sa_tbot_adjust(:,188*4:195*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Temperature (K)','FontSize', fontpt)
       %xticks([24:24:24*18]);
       set(ha(3),'xticklabel',[])

       axes(ha(4));

       ax4=plot(Sa_shum(:,188*4:195*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_shum_adjust(:,188*4:195*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Sa_shum(:,188*4:195*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Sa_shum_adjust(:,188*4:195*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Specific Humidity (kg  kg^{-2})','FontSize', fontpt)
       %xticks([24:24:24*18]);
       linkaxes(ha(1:4), 'x');

       %% Plotting standard deviation staistics to same as Figure(4)
       
       figure(5)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=16;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(6,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       j=plot(std(Faxa_swndr,0,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(std(Faxa_swndr_adjust,0,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM4 Mean', 'Adjusted CAM4 Mean','FontSize',fontpt);
       ylabel('Near-IR, Direct (W  m^{-2})','FontSize', fontpt)
       title(['Standard Deviation: CAM4 scaling against tower NR1 Met: ', yearstr{ii}] , 'FontSize',fontpt);
       set(ha(1),'xticklabel',[])

       axes(ha(2));

       plot(std(Faxa_rain,0,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       plot(std(Faxa_rainl_adjust,0,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Rain (mm  s^{-1})','FontSize', fontpt)
       set(ha(2),'xticklabel',[])

       axes(ha(3));

       plot(std(Faxa_snow,0,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       plot(std(Faxa_snowl_adjust,0,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Snow (mm  s^{-1})','FontSize', fontpt)
       set(ha(3),'xticklabel',[])

       axes(ha(4));

       plot(std(Sa_tbot,0,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       plot(std(Sa_tbot_adjust,0,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Temperature (K)','FontSize', fontpt)
       set(ha(4),'xticklabel',[])

       axes(ha(5));

       plot(std(Sa_shum,0,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       plot(std(Sa_shum_adjust,0,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Specific Humidity (kg  kg^{-2})','FontSize', fontpt)
       set(ha(5),'xticklabel',[])

       axes(ha(6));

       plot(std(Sa_u,0,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       plot(std(Sa_u_adjust,0,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Wind (m s^{-1})','FontSize', fontpt)
       linkaxes(ha(1:6), 'x');

       %% Output standard deviation statistics from above Figure here:


       display(' ')
       display('Display the average 2000-2010 standard deviation for each CAM4 atmospheric forcing')
       display(' ')

       display('Average SW standard deviation for CAM4, adjusted CAM4 and % reduction')
      a= mean(std(Faxa_swndr,0,1)')
      b= mean(std(Faxa_swndr_adjust,0,1)')
      ((a-b)/a)*100  
       display('Average RAIN standard deviation for CAM4, adjusted CAM4 and % reduction')
      a= mean(std(Faxa_rain,0,1)')
      b= mean(std(Faxa_rainl_adjust,0,1)')
      ((a-b)/a)*100
       display('Average SNOW standard deviation for CAM4, adjusted CAM4 and % reduction')
      a= mean(std(Faxa_snow,0,1)')
      b= mean(std(Faxa_snowl_adjust,0,1)')
      ((a-b)/a)*100
       display('Average TEMP standard deviation for CAM4, adjusted CAM4 and % reduction')
      a= mean(std(Sa_tbot,0,1)')
      b= mean(std(Sa_tbot_adjust,0,1)')
      ((a-b)/a)*100
       display('Average SHUM standard deviation for CAM4, adjusted CAM4 and % reduction')
      a= mean(std(Sa_shum,0,1)')
      b= mean(std(Sa_shum_adjust,0,1)')
      ((a-b)/a)*100
       display('Average WIND standard deviation for CAM4, adjusted CAM4 and % reduction')
       a=mean(std(Sa_u,0,1)')
       b=mean(std(Sa_u_adjust,0,1)')
       ((a-b)/a)*100

display('  ')
display('Press enter to proceed to more diagnosic plots..')
display('  ')
pause
close all


            figure(5)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       plot(Faxa_swndr(:,32*4:45*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndr_adjust(:,32*4:45*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       j=plot(mean(Faxa_swndr(:,32*4:45*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       k=plot(mean(Faxa_swndr_adjust(:,32*4:45*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM4 Mean', 'Adjusted CAM4 Mean','FontSize',fontpt);
       ylabel('Near-IR, Direct (W  m^{-2})','FontSize', fontpt)
       title(['CAM4 scaling against tower NR1 Met: February ', yearstr{ii}] , 'FontSize',fontpt);
       xticks([6:6:6*13]);

       axes(ha(2));

       plot(Faxa_snow(:,32*4:45*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_snowl_adjust(:,32*4:45*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_snow(:,32*4:45*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_snowl_adjust(:,32*4:45*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Snow (mm  s^{-1})','FontSize', fontpt)
       xticks([6:6:6*13]);

       axes(ha(3));

       plot(Sa_tbot(:,32*4:45*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_tbot_adjust(:,32*4:45*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Sa_tbot(:,32*4:45*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Sa_tbot_adjust(:,32*4:45*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Temperature (K)','FontSize', fontpt)
       xticks([6:6:6*13]);
       axes(ha(4));

       plot(Sa_shum(:,32*4:45*4)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_shum_adjust(:,32*4:45*4)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Sa_shum(:,32*4:45*4),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Sa_shum_adjust(:,32*4:45*4),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Specific Humidity (kg  kg^{-2})','FontSize', fontpt)
       xticks([6:6:6*13]);

       linkaxes(ha(1:4), 'x');



       figure(6)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(5,1, [.02,.08],[.10,.06],[.08,.06]);

       axes(ha(1));

       j=plot(mean(Faxa_swndr,1)+ ...
              mean(Faxa_swndf,1)+ ...
              mean(Faxa_swvdr,1)+ ...
              mean(Faxa_swvdf,1), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(sw_6hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);
       title(['CAM4 scaling against tower NR1 Met ', yearstr{ii}] , 'FontSize',fontpt);
       ylabel('Total Shortwave (W m^{-2})','FontSize', fontpt)


       axes(ha(2));

       ax1=plot(Faxa_swndr', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndr_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       j=plot(mean(Faxa_swndr,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       k=plot(mean(Faxa_swndr_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM4 Mean', 'Adjusted CAM4 Mean','FontSize',fontpt);
       ylabel('Near-IR, Direct (W m^{-2})','FontSize', fontpt)

       axes(ha(3));

       ax2=plot(Faxa_swndf', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndf_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_swndf,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_swndf_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Near-IR, Diffuse (W m^{-2})','FontSize', fontpt)

       axes(ha(4));

       ax3=plot(Faxa_swvdr', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swvdr_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_swvdr,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_swvdr_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Visible, Direct (W m^{-2})','FontSize', fontpt)

       axes(ha(5));

       ax4=plot(Faxa_swvdf', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swvdf_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_swvdf,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_swvdf_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Visible, Diffuse (W m^{-2})','FontSize', fontpt)
       linkaxes(ha(1:5), 'x');



       figure(7)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);

       axes(ha(1));
       j=plot(mean(Faxa_rainl,1)'+ ...
              mean(Faxa_rainc,1)'+ ...
              mean(Faxa_snowl,1)'+ ...
              mean(Faxa_snowc,1)', ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ppt_6hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);
       title(['CAM4 scaling against tower NR1 Met ', yearstr{ii}] , 'FontSize',fontpt);
       ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)
       ylim([0 2*10^-3])

       axes(ha(2));
       j=plot(mean(Sa_tbot,1), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ta_6hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM4 Ensemble Mean', 'Tower Met','FontSize',fontpt);
       ylabel('Temperature (K)','FontSize', fontpt)

       axes(ha(3));
       plot(Faxa_rain', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_rainl_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_rain,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_rainl_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Rain (mm  s^{-1})','FontSize', fontpt)
       ylim([0 2*10^-3])

       axes(ha(4));

       plot(Faxa_snow', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_snowl_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_snow,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_snowl_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Snow (mm  s^{-1})','FontSize', fontpt)
       ylim([0 2*10^-3])
       linkaxes(ha(1:4), 'x');


       figure(8)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(5,1, [.02,.08],[.10,.06],[.08,.06]);

       axes(ha(1));

       j=plot(mean(Sa_u,1), ...
              '--', 'color',rgb('darksalmon'), 'linewidth', mean_width); hold on
       k=plot(mean(Sa_v,1), ...
              '-', 'color',rgb('darksalmon'), 'linewidth', mean_width);
       l=plot( (mean(Sa_u,1).^2+ ...
               mean(Sa_v,1).^2).^0.5, ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width);

       m=plot(wind_6hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k,l,m], 'Zonal, CAM4 Ens Mean', 'Meridional, CAM4 Ens Mean', ...
                     'Total, CAM4 Ens Mean','Tower Met','FontSize',fontpt);

       ylabel('Wind (m s^{-1})','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr{ii}] , 'FontSize',fontpt);
       axes(ha(2));

       ax1=plot(Sa_u', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_u_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);

       j=plot(mean(Sa_u,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       k=plot(mean(Sa_u_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM4 Mean', 'Adjusted CAM4 Mean','FontSize',fontpt);
       ylabel('Zonal wind (m s^{-1})','FontSize', fontpt)
       axes(ha(3));

       ax2=plot(Sa_v', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_v_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Sa_v,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Sa_v_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Meridional wind (m s^{-1})','FontSize', fontpt)
       axes(ha(4));

       ax3=plot(Sa_pbot', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_pbot_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Sa_pbot,1), '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Sa_pbot_adjust,1), '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Pressure (Pa)','FontSize', fontpt)

       axes(ha(5));

       ax4=plot(Faxa_lwdn', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_lwdn_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Faxa_lwdn,1), '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Faxa_lwdn_adjust,1), '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Longwave Radiation (W m^{-2})','FontSize', fontpt)
       linkaxes(ha(1:5), 'x');




display('  ')
display('Press enter to proceed to CAM4 HEAT MAP plots..')
display('  ')
pause
close all


     % Generate Heat Maps of the CAM4 met forcing relationships
     % Are the relationships fundamentally different?


       figure(6)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=16;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(6,6, [.005,.005],[.06,.06],[.08,.01]);
       
       axes(ha(31));
       % Generate 2 column matrix Faxa_swndr --->(ens,time), (80,8760)
      
       
       SW=Faxa_swndr'+Faxa_swndf'+Faxa_swvdr'+Faxa_swvdf';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       X=[SW_vector,SW_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       %[N,C]=hist3(X,'Nedges',{0:10:1100 0:10:1100})
       [N31,C31]=hist3(X,'Edges',{linspace(0,800,20) linspace(0,800,20)}); hold on
       contourf(C31{1},C31{2},N31)

       colorbar
       caxis([0 10000])
       view(2)
       xlabel('SW (W m^{-2})','FontSize', fontpt)
       ylabel('SW (W m^{-2})','FontSize', fontpt)
      % xlim([0 1100])
      % ylim([0 1100])
       

       axes(ha(32));

       SW=Faxa_swndr'+Faxa_swndf'+Faxa_swvdr'+Faxa_swvdf';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       LW=Faxa_lwdn';
       LW_vector=reshape(LW,numel(LW),1);
       X=[SW_vector,LW_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %Example: {0:10:100 0:50:500}
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N32,C32]=hist3(X,'Edges',{linspace(0,800,20) linspace(125,425,20)}); hold on
       contourf(C32{1},C32{2},N32)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(32),'yticklabel',[])
       xlabel('LW (W m^{-2})','FontSize', fontpt)
       %xlim([125 425])
       %ylim([0 1100])

      axes(ha(33));

       SW=Faxa_swndr'+Faxa_swndf'+Faxa_swvdr'+Faxa_swvdf';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       Pr=Faxa_rainl'+Faxa_rainc'+Faxa_snowl'+Faxa_snowc';
       Pr_vector=reshape(Pr,numel(Pr),1);
       Pr_vector(Pr_vector<0.0000001)=NaN;
       X=[SW_vector,Pr_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N33,C33]=hist3(X,'Edges',{linspace(0,800,20) linspace(0,0.00001,20)}); hold on
       contourf(C33{1},C33{2},N33)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(33),'yticklabel',[])  
       xlabel('Precip (mm s^{-1})','FontSize', fontpt)
       %xlim([0 10*10^(-4)])
       %ylim([0 1100])

       axes(ha(34));

       SW=Faxa_swndr'+Faxa_swndf'+Faxa_swvdr'+Faxa_swvdf';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       Ta=Sa_tbot';
       Ta_vector=reshape(Ta,numel(Ta),1);
       X=[SW_vector,Ta_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N34,C34]=hist3(X,'Edges',{linspace(0,800,20) linspace(250,305,20)}); hold on
       contourf(C34{1},C34{2},N34)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(34),'yticklabel',[])  
       xlabel('Temp (K)','FontSize', fontpt)
       %xlim([250 305])
       %ylim([0 1100])

       axes(ha(35));

       SW=Faxa_swndr'+Faxa_swndf'+Faxa_swvdr'+Faxa_swvdf';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       Shum=Sa_shum;
       Shum_vector=reshape(Shum,numel(Shum),1);
       X=[SW_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N35,C35]=hist3(X,'Edges',{linspace(0,800,20) linspace(0,0.0175,20)}); hold on
       contourf(C35{1},C35{2},N35)
       colorbar
       caxis([0 3000])
       view(2)
       set(ha(35),'yticklabel',[])
       xlabel('Shum (kg kg^{-1})','FontSize', fontpt)
       %xlim([0 17.5*10^(-3)])
       %ylim([0 1100])

       axes(ha(36));

       SW=Faxa_swndr'+Faxa_swndf'+Faxa_swvdr'+Faxa_swvdf';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       Wind=((Sa_u.^2+Sa_v.^2).^0.5)';
       Wind_vector=reshape(Wind,numel(Wind),1);
       X=[SW_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N36,C36]=hist3(X,'Edges',{linspace(0,800,20) linspace(0,20,20)}); hold on
       contourf(C36{1},C36{2},N36)
       colorbar
       caxis([0 3000])
       view(2)
       set(ha(36),'yticklabel',[])
       xlabel('Wind (m s^{-1})','FontSize', fontpt)
       %xlim([0 20])
       %ylim([0 1100])

       % Next row


       axes(ha(26));

       X=[LW_vector,LW_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N26,C26]=hist3(X,'Edges',{linspace(125,425,20) linspace(125,425,20)}); hold on
       contourf(C26{1},C26{2},N26)
       colorbar
       caxis([0 15000])
       view(2)
       set(ha(26),'xticklabel',[])
       ylabel('LW (W m^{-2})','FontSize', fontpt)
       %xlim([125 425])
       %ylim([125 425])

       axes(ha(27));

       X=[LW_vector,Pr_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N27,C27]=hist3(X,'Edges',{linspace(125,425,20) linspace(0,0.00001,20)}); hold on
       contourf(C27{1},C27{2},N27)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(27),'yticklabel',[])
       set(ha(27),'yticklabel',[])
       %xlim([0 10*10^(-4)])
       %ylim([125 425])       

       axes(ha(28));

       X=[LW_vector,Ta_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N28,C28]=hist3(X,'Edges',{linspace(125,425,20) linspace(250,305,20)}); hold on
       contourf(C28{1},C28{2},N28)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(28),'yticklabel',[])
       set(ha(28),'xticklabel',[])
       %xlim([250 305])
       %ylim([125 425])

       axes(ha(29));

       X=[LW_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N29,C29]=hist3(X,'Edges',{linspace(125,425,20) linspace(0,0.0175,20)}); hold on
       contourf(C29{1},C29{2},N29)
       colorbar
       caxis([0 2500])
       view(2)
       set(ha(29),'yticklabel',[])
       set(ha(29),'xticklabel',[])
     
       %xlim([0 17.5*10^(-3)])
       %ylim([125 425])       



       axes(ha(30));

       X=[LW_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N30,C30]=hist3(X,'Edges',{linspace(125,425,20) linspace(0,20,20)}); hold on
       contourf(C30{1},C30{2},N30)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(30),'yticklabel',[])
       set(ha(30),'xticklabel',[])

       %xlim([0 20])
       %ylim([125 425])

       %Next row


       axes(ha(21));

       X=[Pr_vector,Pr_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N21,C21]=hist3(X,'Edges',{linspace(0,0.00001,20) linspace(0,0.00001,20)}); hold on
       contourf(C21{1},C21{2},N21)
       colorbar
       caxis([0 5000])
       view(2)
       ylabel('Precip (mm s^{-1})','FontSize', fontpt)
       set(ha(21),'xticklabel',[])
       %xlim([0 10*10^(-4)])
       %ylim([0 10*10^(-4)])
       
       axes(ha(22));

       X=[Pr_vector,Ta_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N22,C22]=hist3(X,'Edges',{linspace(0,0.00001,20) linspace(250,305,20)}); hold on
       contourf(C22{1},C22{2},N22)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(22),'yticklabel',[])
       set(ha(22),'xticklabel',[])
       %xlim([250 305])
       %ylim([0 10*10^(-4)])


       axes(ha(23));

       X=[Pr_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N23,C23]=hist3(X,'Edges',{linspace(0,0.00001,20) linspace(0,0.0175,20)}); hold on
       contourf(C23{1},C23{2},N23)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(23),'yticklabel',[])
       set(ha(23),'xticklabel',[])
       %xlim([0 17.5*10^(-3)])
       %ylim([0 10*10^(-4)])

       axes(ha(24));

       X=[Pr_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N24,C24]=hist3(X,'Edges',{linspace(0,0.00001,20) linspace(0,20,20)}); hold on 
       contourf(C24{1},C24{2},N24)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(24),'yticklabel',[])
       set(ha(24),'xticklabel',[])
       %xlim([0 20])
       %ylim([0 10*10^(-4)])
    
       % next row
       axes(ha(16));

       X=[Ta_vector,Ta_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N16,C16]=hist3(X,'Edges',{linspace(250,305,20) linspace(250,305,20)}); hold on
       contourf(C16{1},C16{2},N16)
       colorbar
       caxis([0 10000])
       view(2)
       ylabel('Temp (K)','FontSize', fontpt)
       set(ha(16),'xticklabel',[])
       %xlim([250 305])
       %ylim([250 305])

       axes(ha(17));

       X=[Ta_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N17,C17]=hist3(X,'Edges',{linspace(250,305,20) linspace(0,0.0175,20)}); hold on
       contourf(C17{1},C17{2},N17)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(17),'yticklabel',[])
       set(ha(17),'xticklabel',[])
       %ylim([250 305])
       %xlim([0 17.5*10^(-3)])

       axes(ha(18));

       X=[Ta_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N18,C18]=hist3(X,'Edges',{linspace(250,305,20) linspace(0,20,20)}); hold on
       contourf(C18{1},C18{2},N18)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(18),'yticklabel',[])
       set(ha(18),'xticklabel',[])
       %ylim([250 305])
       %xlim([0 20])

       axes(ha(11));

       X=[Shum_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N11,C11]=hist3(X,'Edges',{linspace(0,0.0175,20) linspace(0,0.0175,20)}); hold on
       contourf(C11{1},C11{2},N11)
       colorbar
       caxis([0 15000])
       view(2)
       ylabel('Shum (kg kg^{-1})','FontSize', fontpt)
       set(ha(11),'xticklabel',[])
       %xlim([0 17.5*10^(-3)])
       %ylim([0 17.5*10^(-3)])

       axes(ha(12));

       X=[Shum_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N12,C12]=hist3(X,'Edges',{linspace(0,0.0175,20) linspace(0,20,20)}); hold on
       contourf(C12{1},C12{2},N12)
       colorbar
       caxis([0 2500])
       view(2)
       set(ha(12),'yticklabel',[])
       set(ha(12),'xticklabel',[])
       %xlim([0 20])
       %ylim([0 17.5*10^(-3)])

       axes(ha(6));

       X=[Wind_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N6,C6]=hist3(X,'Edges',{linspace(0,20,20) linspace(0,20,20)}); hold on
       contourf(C6{1},C6{2},N6)
       colorbar
       caxis([0 10000])
       view(2)
       set(ha(6),'xticklabel',[])
       ylabel('Wind (m s^{-1})','FontSize', fontpt)
       %xlim([0 20])
       %ylim([0 20])




%% Switch to the CAM4 bias adjusted  correction structure


       figure(7)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=16;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(6,6, [.005,.005],[.06,.06],[.08,.01]);
 
       axes(ha(31));
       SW=Faxa_swndr_adjust'+Faxa_swndf_adjust'+Faxa_swvdr_adjust'+Faxa_swvdf_adjust';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       X=[SW_vector,SW_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       %[N,C]=hist3(X,'Nedges',{0:10:1100 0:10:1100})
       [N31_adj,C31_adj]=hist3(X,'Edges',{linspace(0,800,20) linspace(0,800,20)}); hold on
       contourf(C31_adj{1},C31_adj{2},N31_adj)

       colorbar
       caxis([0 10000])
       view(2)
       xlabel('SW adjust (W m^{-2})','FontSize', fontpt)
       ylabel('SW adjust (W m^{-2})','FontSize', fontpt)
      % xlim([0 1100])
      % ylim([0 1100])


       axes(ha(32));

       SW=Faxa_swndr_adjust'+Faxa_swndf_adjust'+Faxa_swvdr_adjust'+Faxa_swvdf_adjust';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       LW=Faxa_lwdn_adjust';
       LW_vector=reshape(LW,numel(LW),1);
       X=[SW_vector,LW_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %Example: {0:10:100 0:50:500}
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N32_adj,C32_adj]=hist3(X,'Edges',{linspace(0,800,20) linspace(125,425,20)}); hold on
       contourf(C32_adj{1},C32_adj{2},N32_adj)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(32),'yticklabel',[])
       xlabel('LW adjust (W m^{-2})','FontSize', fontpt)
       %xlim([125 425])
       %ylim([0 1100])

      axes(ha(33));

       SW=Faxa_swndr_adjust'+Faxa_swndf_adjust'+Faxa_swvdr_adjust'+Faxa_swvdf_adjust';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector==0)=NaN;
       SW_vector(SW_vector<5)=NaN;
       Pr=Faxa_rainl_adjust'+Faxa_rainc_adjust'+Faxa_snowl_adjust'+Faxa_snowc_adjust';
       Pr_vector=reshape(Pr,numel(Pr),1);
       Pr_vector(Pr_vector<0.0000001)=NaN;
       X=[SW_vector,Pr_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N33_adj,C33_adj]=hist3(X,'Edges',{linspace(0,800,20) linspace(0,0.00001,20)}); hold on
       contourf(C33_adj{1},C33_adj{2},N33_adj)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(33),'yticklabel',[])
       xlabel('Precip adjust (mm s^{-1})','FontSize', fontpt)
       %xlim([0 10*10^(-4)])
       %ylim([0 1100])

       axes(ha(34));

       SW=Faxa_swndr_adjust'+Faxa_swndf_adjust'+Faxa_swvdr_adjust'+Faxa_swvdf_adjust';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       Ta=Sa_tbot_adjust';
       Ta_vector=reshape(Ta,numel(Ta),1);
       X=[SW_vector,Ta_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N34_adj,C34_adj]=hist3(X,'Edges',{linspace(0,800,20) linspace(250,305,20)}); hold on
       contourf(C34_adj{1},C34_adj{2},N34_adj)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(34),'yticklabel',[])
       xlabel('Temp adjust (K)','FontSize', fontpt)
       %xlim([250 305])
       %ylim([0 1100])

       axes(ha(35));

       SW=Faxa_swndr_adjust'+Faxa_swndf_adjust'+Faxa_swvdr_adjust'+Faxa_swvdf_adjust';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       Shum=Sa_shum_adjust';
       Shum_vector=reshape(Shum,numel(Shum),1);
       X=[SW_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N35_adj,C35_adj]=hist3(X,'Edges',{linspace(0,800,20) linspace(0,0.0175,20)}); hold on
       contourf(C35_adj{1},C35_adj{2},N35_adj)
       colorbar
       caxis([0 3000])
       view(2)
       set(ha(35),'yticklabel',[])
       xlabel('Shum adjust (kg kg^{-1})','FontSize', fontpt)
       %xlim([0 17.5*10^(-3)])
       %ylim([0 1100])

       axes(ha(36));

       SW=Faxa_swndr_adjust'+Faxa_swndf_adjust'+Faxa_swvdr_adjust'+Faxa_swvdf_adjust';
       SW_vector=reshape(SW,numel(SW),1);
       SW_vector(SW_vector<5)=NaN;
       Wind=((Sa_u_adjust.^2+Sa_v_adjust.^2).^0.5)';
       Wind_vector=reshape(Wind,numel(Wind),1);
       X=[SW_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N36_adj,C36_adj]=hist3(X,'Edges',{linspace(0,800,20) linspace(0,20,20)}); hold on
       contourf(C36_adj{1},C36_adj{2},N36_adj)
       colorbar
       caxis([0 3000])
       view(2)
       set(ha(36),'yticklabel',[])
       xlabel('Wind (m s^{-1})','FontSize', fontpt)
       %xlim([0 20])
       %ylim([0 1100])

       % Next row


       axes(ha(26));

       X=[LW_vector,LW_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N26_adj,C26_adj]=hist3(X,'Edges',{linspace(125,425,20) linspace(125,425,20)}); hold on
       contourf(C26_adj{1},C26_adj{2},N26_adj)
       colorbar
       caxis([0 15000])
       view(2)
       set(ha(26),'xticklabel',[])
       ylabel('LW adjust (W m^{-2})','FontSize', fontpt)
       %xlim([125 425])
       %ylim([125 425])

       axes(ha(27));

       X=[LW_vector,Pr_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N27_adj,C27_adj]=hist3(X,'Edges',{linspace(125,425,20) linspace(0,0.00001,20)}); hold on
       contourf(C27_adj{1},C27_adj{2},N27_adj)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(27),'yticklabel',[])
       set(ha(27),'yticklabel',[])
       %xlim([0 10*10^(-4)])
       %ylim([125 425])       

       axes(ha(28));

       X=[LW_vector,Ta_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N28_adj,C28_adj]=hist3(X,'Edges',{linspace(125,425,20) linspace(250,305,20)}); hold on
       contourf(C28_adj{1},C28_adj{2},N28_adj)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(28),'yticklabel',[])
       set(ha(28),'xticklabel',[])
       %xlim([250 305])
       %ylim([125 425])

       axes(ha(29));

       X=[LW_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N29_adj,C29_adj]=hist3(X,'Edges',{linspace(125,425,20) linspace(0,0.0175,20)}); hold on
       contourf(C29_adj{1},C29_adj{2},N29_adj)
       colorbar
       caxis([0 2500])
       view(2)
       set(ha(29),'yticklabel',[])
       set(ha(29),'xticklabel',[])

       %xlim([0 17.5*10^(-3)])
       %ylim([125 425])       



       axes(ha(30));

       X=[LW_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N30_adj,C30_adj]=hist3(X,'Edges',{linspace(125,425,20) linspace(0,20,20)}); hold on
       contourf(C30_adj{1},C30_adj{2},N30_adj)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(30),'yticklabel',[])
       set(ha(30),'xticklabel',[])

       %xlim([0 20])
       %ylim([125 425])

       %Next row


       axes(ha(21));

       X=[Pr_vector,Pr_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N21_adj,C21_adj]=hist3(X,'Edges',{linspace(0,0.00001,20) linspace(0,0.00001,20)}); hold on
       contourf(C21_adj{1},C21_adj{2},N21_adj)
       colorbar
       caxis([0 5000])
       view(2)
       ylabel('Precip adjust (mm s^{-1})','FontSize', fontpt)
       set(ha(21),'xticklabel',[])
       %xlim([0 10*10^(-4)])
       %ylim([0 10*10^(-4)])

       axes(ha(22));

       X=[Pr_vector,Ta_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N22_adj,C22_adj]=hist3(X,'Edges',{linspace(0,0.00001,20) linspace(250,305,20)}); hold on
       contourf(C22_adj{1},C22_adj{2},N22_adj)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(22),'yticklabel',[])
       set(ha(22),'xticklabel',[])
       %xlim([250 305])
       %ylim([0 10*10^(-4)])


       axes(ha(23));

       X=[Pr_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N23_adj,C23_adj]=hist3(X,'Edges',{linspace(0,0.00001,20) linspace(0,0.0175,20)}); hold on
       contourf(C23_adj{1},C23_adj{2},N23_adj)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(23),'yticklabel',[])
       set(ha(23),'xticklabel',[])
       %xlim([0 17.5*10^(-3)])
       %ylim([0 10*10^(-4)])

       axes(ha(24));

       X=[Pr_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N24_adj,C24_adj]=hist3(X,'Edges',{linspace(0,0.00001,20) linspace(0,20,20)}); hold on
       contourf(C24_adj{1},C24_adj{2},N24_adj)
       colorbar
       caxis([0 1000])
       view(2)
       set(ha(24),'yticklabel',[])
       set(ha(24),'xticklabel',[])
       %xlim([0 20])
       %ylim([0 10*10^(-4)])

       % next row
       axes(ha(16));

       X=[Ta_vector,Ta_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N16_adj,C16_adj]=hist3(X,'Edges',{linspace(250,305,20) linspace(250,305,20)}); hold on
       contourf(C16_adj{1},C16_adj{2},N16_adj)
       colorbar
       caxis([0 10000])
       view(2)
       ylabel('Temp adjust (K)','FontSize', fontpt)
       set(ha(16),'xticklabel',[])
       %xlim([250 305])
       %ylim([250 305])

       axes(ha(17));

       X=[Ta_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N17_adj,C17_adj]=hist3(X,'Edges',{linspace(250,305,20) linspace(0,0.0175,20)}); hold on
       contourf(C17_adj{1},C17_adj{2},N17_adj)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(17),'yticklabel',[])
       set(ha(17),'xticklabel',[])
       %ylim([250 305])
       %xlim([0 17.5*10^(-3)])

       axes(ha(18));

       X=[Ta_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N18_adj,C18_adj]=hist3(X,'Edges',{linspace(250,305,20) linspace(0,20,20)}); hold on
       contourf(C18_adj{1},C18_adj{2},N18_adj)
       colorbar
       caxis([0 2000])
       view(2)
       set(ha(18),'yticklabel',[])
       set(ha(18),'xticklabel',[])
       %ylim([250 305])
       %xlim([0 20])

       axes(ha(11));

       X=[Shum_vector,Shum_vector];
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %hist3(X,'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N11_adj,C11_adj]=hist3(X,'Edges',{linspace(0,0.0175,20) linspace(0,0.0175,20)}); hold on
       contourf(C11_adj{1},C11_adj{2},N11_adj)
       colorbar
       caxis([0 15000])
       view(2)
       ylabel('Shum adjust (kg kg^{-1})','FontSize', fontpt)
       set(ha(11),'xticklabel',[])
       %xlim([0 17.5*10^(-3)])
       %ylim([0 17.5*10^(-3)])


       axes(ha(12));

       X=[Shum_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N12_adj,C12_adj]=hist3(X,'Edges',{linspace(0,0.0175,20) linspace(0,20,20)}); hold on
       contourf(C12_adj{1},C12_adj{2},N12_adj)
       colorbar
       caxis([0 2500])
       view(2)
       set(ha(12),'yticklabel',[])
       set(ha(12),'xticklabel',[])
       %xlim([0 20])
       %ylim([0 17.5*10^(-3)])

       axes(ha(6));

       X=[Wind_vector,Wind_vector];
       %hist3(X,'Cdatamode','auto'); hold on
       %hist3(X,'Nbins',[20 20],'Cdatamode','auto'); hold on
       %[N,C]=hist3(X,'Nbins',[20 20])
       [N6_adj,C6_adj]=hist3(X,'Edges',{linspace(0,20,20) linspace(0,20,20)}); hold on
       contourf(C6_adj{1},C6_adj{2},N6_adj)
       colorbar
       caxis([0 10000])
       view(2)
       set(ha(6),'xticklabel',[])
       ylabel('Wind adjust (m s^{-1})','FontSize', fontpt)
       %xlim([0 20])
       %ylim([0 20])

      %% Calculate the grid count differences for the original and adjusted CAM4

     N6_diff=N6_adj-N6; C6_diff1=C6_adj{1}-C6{1}; C6_diff2=C6_adj{2}-C6{2};

     N12_diff=N12_adj-N12; C12_diff1=C12_adj{1}-C12{1}; C12_diff2=C12_adj{2}-C12{2};
     N11_diff=N11_adj-N6; C11_diff1=C11_adj{1}-C11{1}; C11_diff2=C11_adj{2}-C11{2};

     N16_diff=N16_adj-N16; C16_diff1=C16_adj{1}-C16{1}; C16_diff2=C16_adj{2}-C16{2};
     N17_diff=N17_adj-N17; C17_diff1=C17_adj{1}-C17{1}; C17_diff2=C17_adj{2}-C17{2};
     N18_diff=N18_adj-N18; C18_diff1=C18_adj{1}-C18{1}; C18_diff2=C18_adj{2}-C18{2};

     N21_diff=N21_adj-N21; C21_diff1=C21_adj{1}-C21{1}; C21_diff2=C21_adj{2}-C21{2};
     N22_diff=N22_adj-N22; C22_diff1=C22_adj{1}-C22{1}; C22_diff2=C22_adj{2}-C22{2};
     N23_diff=N23_adj-N23; C23_diff1=C23_adj{1}-C23{1}; C23_diff2=C23_adj{2}-C23{2};
     N24_diff=N24_adj-N24; C24_diff1=C24_adj{1}-C24{1}; C24_diff2=C24_adj{2}-C24{2};

     N26_diff=N26_adj-N26; C26_diff1=C26_adj{1}-C26{1}; C26_diff2=C26_adj{2}-C26{2};
     N27_diff=N27_adj-N27; C27_diff1=C27_adj{1}-C27{1}; C27_diff2=C27_adj{2}-C27{2};
     N28_diff=N28_adj-N28; C28_diff1=C28_adj{1}-C28{1}; C28_diff2=C28_adj{2}-C28{2};
     N29_diff=N29_adj-N29; C29_diff1=C29_adj{1}-C29{1}; C29_diff2=C29_adj{2}-C29{2};
     N30_diff=N30_adj-N30; C30_diff1=C30_adj{1}-C30{1}; C30_diff2=C30_adj{2}-C30{2};
   
     N31_diff=N31_adj-N31; C31_diff1=C31_adj{1}-C31{1}; C31_diff2=C31_adj{2}-C31{2};
     N32_diff=N32_adj-N32; C32_diff1=C32_adj{1}-C32{1}; C32_diff2=C32_adj{2}-C32{2};
     N33_diff=N33_adj-N33; C33_diff1=C33_adj{1}-C33{1}; C33_diff2=C33_adj{2}-C33{2};
     N34_diff=N34_adj-N34; C34_diff1=C34_adj{1}-C34{1}; C34_diff2=C34_adj{2}-C34{2};
     N35_diff=N35_adj-N35; C35_diff1=C35_adj{1}-C35{1}; C35_diff2=C35_adj{2}-C35{2};
     N36_diff=N36_adj-N36; C36_diff1=C36_adj{1}-C36{1}; C36_diff2=C36_adj{2}-C36{2};


       figure(8)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=16;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(6,6, [.005,.005],[.06,.06],[.08,.01]);

       axes(ha(31));
       contourf(C31{1},C31{2},N31_diff)

       
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 10000])
       view(2)
       xlabel('SW Diff (W m^{-2})','FontSize', fontpt)
       ylabel('SW Diff (W m^{-2})','FontSize', fontpt)
      % xlim([0 1100])
      % ylim([0 1100])


       axes(ha(32));

       contourf(C32{1},C32{2},N32_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 3000])
       view(2)
       set(ha(32),'yticklabel',[])
       xlabel('LW Diff (W m^{-2})','FontSize', fontpt)
       %xlim([125 425])
       %ylim([0 1100])

      axes(ha(33));

       contourf(C33{1},C33{2},N33_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 3000])
       view(2)
       set(ha(33),'yticklabel',[])
       xlabel('Precip adjust (mm s^{-1})','FontSize', fontpt)
       %xlim([0 10*10^(-4)])
       %ylim([0 1100])

       axes(ha(34));

       contourf(C34{1},C34{2},N34_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 3000])
       view(2)
       set(ha(34),'yticklabel',[])
       xlabel('Temp Diff (K)','FontSize', fontpt)
       %xlim([250 305])
       %ylim([0 1100])

       axes(ha(35));

       contourf(C35{1},C35{2},N35_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 5000])
       view(2)
       set(ha(35),'yticklabel',[])
       xlabel('Shum Diff (kg kg^{-1})','FontSize', fontpt)
       %xlim([0 17.5*10^(-3)])
       %ylim([0 1100])

       axes(ha(36));

       contourf(C36{1},C36{2},N36_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 5000])
       view(2)
       set(ha(36),'yticklabel',[])
       xlabel('Wind Diff (m s^{-1})','FontSize', fontpt)
       %xlim([0 20])
       %ylim([0 1100])

       % Next row


       axes(ha(26));

       contourf(C26{1},C26{2},N26_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 20000])
       view(2)
       set(ha(26),'xticklabel',[])
       ylabel('LW Diff (W m^{-2})','FontSize', fontpt)
       %xlim([125 425])
       %ylim([125 425])

       axes(ha(27));

       contourf(C27{1},C27{2},N27_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 7000])
       view(2)
       set(ha(27),'yticklabel',[])
       set(ha(27),'yticklabel',[])
       %xlim([0 10*10^(-4)])
       %ylim([125 425])       

       axes(ha(28));

       contourf(C28{1},C28{2},N28_diff)
       colormap(b2r(-10000,10000));
       colorbar
       %caxis([0 20000])
       view(2)
       set(ha(28),'yticklabel',[])
       set(ha(28),'xticklabel',[])
       %xlim([250 305])
       %ylim([125 425])

       axes(ha(29));

       contourf(C29{1},C29{2},N29_diff)
       colormap(b2r(-10000,10000));
       colorbar
       %caxis([0 20000])
       view(2)
       set(ha(29),'yticklabel',[])
       set(ha(29),'xticklabel',[])

       %xlim([0 17.5*10^(-3)])
       %ylim([125 425])       



       axes(ha(30));

       contourf(C30{1},C30{2},N30_diff)
       colormap(b2r(-10000,10000));
       colorbar
       %caxis([0 12000])
       view(2)
       set(ha(30),'yticklabel',[])
       set(ha(30),'xticklabel',[])

       %xlim([0 20])
       %ylim([125 425])

       %Next row


       axes(ha(21));

       contourf(C21{1},C21{2},N21_diff)
       colormap(b2r(-5000,5000)); 
       colorbar
       %caxis([0 20000])
       view(2)
       ylabel('Precip Diff (mm s^{-1})','FontSize', fontpt)
       set(ha(21),'xticklabel',[])
       %xlim([0 10*10^(-4)])
       %ylim([0 10*10^(-4)])

       axes(ha(22));

       contourf(C22{1},C22{2},N22_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 5000])
       view(2)
       set(ha(22),'yticklabel',[])
       set(ha(22),'xticklabel',[])
       %xlim([250 305])
       %ylim([0 10*10^(-4)])


       axes(ha(23));

       contourf(C23{1},C23{2},N23_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 5000])
       view(2)
       set(ha(23),'yticklabel',[])
       set(ha(23),'xticklabel',[])
       %xlim([0 17.5*10^(-3)])
       %ylim([0 10*10^(-4)])

       axes(ha(24));

       contourf(C24{1},C24{2},N24_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 5000])
       view(2)
       set(ha(24),'yticklabel',[])
       set(ha(24),'xticklabel',[])
       %xlim([0 20])
       %ylim([0 10*10^(-4)])

       % next row
       axes(ha(16));

       contourf(C16{1},C16{2},N16_diff)
       colormap(b2r(-10000,10000));
       colorbar
       %caxis([0 25000])
       view(2)
       ylabel('Temp Diff (K)','FontSize', fontpt)
       set(ha(16),'xticklabel',[])
       %xlim([250 305])
       %ylim([250 305])

       axes(ha(17));

       contourf(C17{1},C17{2},N17_diff)
       colormap(b2r(-10000,10000));
       colorbar
       %caxis([0 15000])
       view(2)
       set(ha(17),'yticklabel',[])
       set(ha(17),'xticklabel',[])
       %ylim([250 305])
       %xlim([0 17.5*10^(-3)])

       axes(ha(18));

       contourf(C18{1},C18{2},N18_diff)
       colormap(b2r(-10000,10000));
       colorbar
       %caxis([0 12000])
       view(2)
       set(ha(18),'yticklabel',[])
       set(ha(18),'xticklabel',[])
       %ylim([250 305])
       %xlim([0 20])

       axes(ha(11));

       contourf(C11{1},C11{2},N11_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 25000])
       view(2)
       ylabel('Shum Diff (kg kg^{-1})','FontSize', fontpt)
       set(ha(11),'xticklabel',[])
       %xlim([0 17.5*10^(-3)])
       %ylim([0 17.5*10^(-3)])


       axes(ha(12));

       contourf(C12{1},C12{2},N12_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 17000])
       view(2)
       set(ha(12),'yticklabel',[])
       set(ha(12),'xticklabel',[])
       %xlim([0 20])
       %ylim([0 17.5*10^(-3)])

       axes(ha(6));

       contourf(C6{1},C6{2},N6_diff)
       colormap(b2r(-5000,5000));
       colorbar
       %caxis([0 25000])
       view(2)
       set(ha(6),'xticklabel',[])
       ylabel('Wind Diff (m s^{-1})','FontSize', fontpt)
       %xlim([0 20])
       %ylim([0 20])


display('   ')
display(['Press enter to write ' yearstr{ii} ' CAM4 adjusted files to netcdf'])
display('   ')
%pause
close all

stop
end % Diagnostics Loop
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      % Assign the adjusted CAM variables to netcdf files for each ensemble member/ year

          time_var= ncread([path_CAM_forcing enstr{1} '/CAM4_NR1.cpl_' enstr{1} '.ha2x1dx6h.' yearstr{ii} '.nc'],'time');
          dim_time=length(time_var);        % ha2x1dx6h
          time_bnds = ncread([path_CAM_forcing enstr{1} '/CAM4_NR1.cpl_' enstr{1} '.ha2x1dx6h.' yearstr{ii} '.nc'],'time_bnds');
          path_scaled_CAM_forcing = '/glade/work/bmraczka/CAM4_NR1/scaled/';

       
      % Allocate met variables in dimension necessary for assignment to netcdf

        swndr_assign=ones(1,1,length(Faxa_swndr(1,:)))*NaN;  swvdr_assign=ones(1,1,length(Faxa_swvdr(1,:)))*NaN;
        swndf_assign=ones(1,1,length(Faxa_swndf(1,:)))*NaN;  swvdf_assign=ones(1,1,length(Faxa_swvdf(1,:)))*NaN;
        rainc_assign=ones(1,1,length(Faxa_rainc(1,:)))*NaN;  rainl_assign=ones(1,1,length(Faxa_rainl(1,:)))*NaN;
        snowc_assign=ones(1,1,length(Faxa_snowc(1,:)))*NaN; snowl_assign=ones(1,1,length(Faxa_snowl(1,:)))*NaN;
        u_assign=ones(1,1,length(Sa_u(1,:)))*NaN;            v_assign=ones(1,1,length(Sa_v(1,:)))*NaN;
        tbot_assign=ones(1,1,length(Sa_tbot(1,:)))*NaN;      shum_assign=ones(1,1,length(Sa_shum(1,:)))*NaN;
        pbot_assign=ones(1,1,length(Sa_pbot(1,:)))*NaN;      lwdn_assign=ones(1,1,length(Faxa_lwdn(1,:)))*NaN;
        z_assign=ones(1,1,length(Sa_z(1,:)))*NaN;            ptem_assign=ones(1,1,length(Sa_ptem(1,:)))*NaN;
        dens_assign=ones(1,1,length(Sa_dens(1,:)))*NaN;      pslv_assign=ones(1,1,length(Sa_pslv(1,:)))*NaN;

       for jj = 1:80;  %% ensemble loop         
       %for jj = 1:1; %TEST 
        % Assign the met variables to appropriate 3-dimension format
        
        swndr_assign(1,1,:)=Faxa_swndr_adjust(jj,:);                swvdr_assign(1,1,:)=Faxa_swvdr_adjust(jj,:);
        swndf_assign(1,1,:)=Faxa_swndf_adjust(jj,:);                swvdf_assign(1,1,:)=Faxa_swvdf_adjust(jj,:);
        rainc_assign(1,1,:)=Faxa_rainc_adjust(jj,:);                rainl_assign(1,1,:)=Faxa_rainl_adjust(jj,:);
        snowc_assign(1,1,:)=Faxa_snowc_adjust(jj,:);                snowl_assign(1,1,:)=Faxa_snowl_adjust(jj,:);
        u_assign(1,1,:)=Sa_u_adjust(jj,:);                          v_assign(1,1,:)=Sa_v_adjust(jj,:);
        tbot_assign(1,1,:)=Sa_tbot_adjust(jj,:);                    shum_assign(1,1,:)=Sa_shum_adjust(jj,:);
        pbot_assign(1,1,:)=Sa_pbot_adjust(jj,:);                    lwdn_assign(1,1,:)=Faxa_lwdn_adjust(jj,:);
        z_assign(1,1,:)=Sa_z(jj,:);                                 ptem_assign(1,1,:)=Sa_ptem(jj,:);
        dens_assign(1,1,:)=Sa_dens(jj,:);                           pslv_assign(1,1,:)=Sa_pslv(jj,:);



        ncname=[path_scaled_CAM_forcing 'CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'];
        ncid=netcdf.create(ncname,'CLOBBER');
        netcdf.close(ncid);

        nccreate(ncname,'time','Dimensions',{'time',Inf},'Format','classic');
        ncwrite(ncname,'time',time_var)

        ncwriteatt(ncname,'time','units',['days since ' yearstr{ii} '-01-01 00:00:00'])

        ncwriteatt(ncname,'time','calendar','noleap')
        ncwriteatt(ncname,'time','bounds','time_bnds')

        nccreate(ncname,'time_bnds','Dimensions',{'ntb',2,'time',Inf},'Format','classic');
        ncwrite(ncname,'time_bnds',time_bnds)


        nccreate(ncname,'doma_area','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_area',SITE_doma_area)
        ncwriteatt(ncname,'doma_area','_FillValue', 0.)
        ncwriteatt(ncname,'doma_area','units','m^2')
        ncwriteatt(ncname,'doma_area','long_name','undefined')
        ncwriteatt(ncname,'doma_area','standard_name','cell area')
        ncwriteatt(ncname,'doma_area','internal_dname','doma')
        
        nccreate(ncname,'doma_lat','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lat',SITE_lat)
        ncwriteatt(ncname,'doma_lat','_FillValue', 0.)
        ncwriteatt(ncname,'doma_lat','units','degrees north')
        ncwriteatt(ncname,'doma_lat','long_name','undefined')
        ncwriteatt(ncname,'doma_lat','standard_name','latitude')
        ncwriteatt(ncname,'doma_lat','internal_dname','doma')
        
        nccreate(ncname,'doma_lon','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lon',SITE_lon)
        ncwriteatt(ncname,'doma_lon','_FillValue', 0.)
        ncwriteatt(ncname,'doma_lon','units','degrees east')
        ncwriteatt(ncname,'doma_lon','long_name','undefined')
        ncwriteatt(ncname,'doma_lon','standard_name','longitude')
        ncwriteatt(ncname,'doma_lon','internal_dname','doma')

        nccreate(ncname,'doma_mask','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_mask',SITE_doma_mask)
        ncwriteatt(ncname,'doma_mask','_FillValue', 0.)
        ncwriteatt(ncname,'doma_mask','units','unitless')
        ncwriteatt(ncname,'doma_mask','long_name','undefined')
        ncwriteatt(ncname,'doma_mask','standard_name','mask')
        ncwriteatt(ncname,'doma_mask','internal_dname','doma')
        
        nccreate(ncname,'a2x6h_Faxa_swndf','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_swndf',swndf_assign)
        ncwriteatt(ncname,'a2x6h_Faxa_swndf','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Faxa_swndf','units','W m-2')
        ncwriteatt(ncname,'a2x6h_Faxa_swndf','long_name','Diffuse near-infrared incident solar radiation')
        ncwriteatt(ncname,'a2x6h_Faxa_swndf','standard_name','surface_downward_diffuse_shortwave_flux_due_to_near_infrared_radiation')
        ncwriteatt(ncname,'a2x6h_Faxa_swndf','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Faxa_swndf','cell_methods', 'time: mean')

        nccreate(ncname,'a2x6h_Faxa_swndr','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_swndr',swndr_assign)
        ncwriteatt(ncname,'a2x6h_Faxa_swndr','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Faxa_swndr','units','W m-2')
        ncwriteatt(ncname,'a2x6h_Faxa_swndr','long_name','Direct near-infrared incident solar radiation')
        ncwriteatt(ncname,'a2x6h_Faxa_swndr','standard_name','surface_downward_direct_shortwave_flux_due_to_near_infrared_radiation')
        ncwriteatt(ncname,'a2x6h_Faxa_swndr','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Faxa_swndr','cell_methods', 'time: mean') 

        nccreate(ncname,'a2x6h_Faxa_swvdf','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_swvdf',swvdf_assign)
        ncwriteatt(ncname,'a2x6h_Faxa_swvdf','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Faxa_swvdf','units','W m-2')
        ncwriteatt(ncname,'a2x6h_Faxa_swvdf','long_name','Diffuse visible incident solar radiation')
        ncwriteatt(ncname,'a2x6h_Faxa_swvdf','standard_name','surface_downward_diffuse_shortwave_flux_due_to_visible_radiation')
        ncwriteatt(ncname,'a2x6h_Faxa_swvdf','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Faxa_swvdf','cell_methods', 'time: mean')

        nccreate(ncname,'a2x6h_Faxa_swvdr','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_swvdr',swvdr_assign)
        ncwriteatt(ncname,'a2x6h_Faxa_swvdr','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Faxa_swvdr','units','W m-2')
        ncwriteatt(ncname,'a2x6h_Faxa_swvdr','long_name','Direct visible incident solar radiation')
        ncwriteatt(ncname,'a2x6h_Faxa_swvdr','standard_name','surface_downward_direct_shortwave_flux_due_to_visible_radiation')
        ncwriteatt(ncname,'a2x6h_Faxa_swvdr','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Faxa_swvdr','cell_methods', 'time: mean')
        

                                                                                                                
        nccreate(ncname,'a2x6h_Sa_z','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_z',z_assign)
        ncwriteatt(ncname,'a2x6h_Sa_z','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Sa_z','units','m')
        ncwriteatt(ncname,'a2x6h_Sa_z','long_name','Height at the lowest model level')
        ncwriteatt(ncname,'a2x6h_Sa_z','standard_name','height')
        ncwriteatt(ncname,'a2x6h_Sa_z','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Sa_z','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Sa_tbot','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_tbot',tbot_assign)
        ncwriteatt(ncname,'a2x6h_Sa_tbot','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Sa_tbot','units','K')
        ncwriteatt(ncname,'a2x6h_Sa_tbot','long_name','Temperature at the lowest model level')
        ncwriteatt(ncname,'a2x6h_Sa_tbot','standard_name','air_temperature')
        ncwriteatt(ncname,'a2x6h_Sa_tbot','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Sa_tbot','cell_methods','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_ptem','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_ptem',ptem_assign)
        ncwriteatt(ncname,'a2x6h_Sa_ptem','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Sa_ptem','units','K')
        ncwriteatt(ncname,'a2x6h_Sa_ptem','long_name','Potential temperature at the lowest model level')
        ncwriteatt(ncname,'a2x6h_Sa_ptem','standard_name','air_potential_temperature')
        ncwriteatt(ncname,'a2x6h_Sa_ptem','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Sa_ptem','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Sa_shum','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_shum',shum_assign)
        ncwriteatt(ncname,'a2x6h_Sa_shum','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Sa_shum','units','kg kg-1')
        ncwriteatt(ncname,'a2x6h_Sa_shum','long_name','Specific humidity at the lowest model level')
        ncwriteatt(ncname,'a2x6h_Sa_shum','standard_name','specific humidity')
        ncwriteatt(ncname,'a2x6h_Sa_shum','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Sa_shum','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Sa_dens','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_dens',dens_assign)
        ncwriteatt(ncname,'a2x6h_Sa_dens','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Sa_dens','units','kg m-3')
        ncwriteatt(ncname,'a2x6h_Sa_dens','long_name','Air density at the lowest model level')
        ncwriteatt(ncname,'a2x6h_Sa_dens','standard_name','air_density')
        ncwriteatt(ncname,'a2x6h_Sa_dens','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Sa_dens','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Sa_pbot','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_pbot',pbot_assign)
        ncwriteatt(ncname,'a2x6h_Sa_pbot','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Sa_pbot','units','Pa')
        ncwriteatt(ncname,'a2x6h_Sa_pbot','long_name','Pressure at the lowest model level')
        ncwriteatt(ncname,'a2x6h_Sa_pbot','standard_name','air_pressure')
        ncwriteatt(ncname,'a2x6h_Sa_pbot','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Sa_pbot','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Sa_pslv','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_pslv',pslv_assign)
        ncwriteatt(ncname,'a2x6h_Sa_pslv','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Sa_pslv','units','Pa')
        ncwriteatt(ncname,'a2x6h_Sa_pslv','long_name','Sea level pressure')
        ncwriteatt(ncname,'a2x6h_Sa_pslv','standard_name','air_pressure_at_sea_level')
        ncwriteatt(ncname,'a2x6h_Sa_pslv','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Sa_pslv','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Faxa_lwdn','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_lwdn',lwdn_assign)
        ncwriteatt(ncname,'a2x6h_Faxa_lwdn','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Faxa_lwdn','units','W m-2')
        ncwriteatt(ncname,'a2x6h_Faxa_lwdn','long_name','Downward longwave heat flux')
        ncwriteatt(ncname,'a2x6h_Faxa_lwdn','standard_name','downwelling_longwave_flux')
        ncwriteatt(ncname,'a2x6h_Faxa_lwdn','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Faxa_lwdn','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Faxa_rainc','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_rainc',rainc_assign)
        ncwriteatt(ncname,'a2x6h_Faxa_rainc','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Faxa_rainc','units','kg m-2 s-1')
        ncwriteatt(ncname,'a2x6h_Faxa_rainc','long_name','Convective precipitation rate')
        ncwriteatt(ncname,'a2x6h_Faxa_rainc','standard_name','convective_precipitation_flux')
        ncwriteatt(ncname,'a2x6h_Faxa_rainc','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Faxa_rainc','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Faxa_rainl','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_rainl',rainl_assign)
        ncwriteatt(ncname,'a2x6h_Faxa_rainl','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Faxa_rainl','units','kg m-2 s-1')
        ncwriteatt(ncname,'a2x6h_Faxa_rainl','long_name','Large-scale (stable) precipitation rate')
        ncwriteatt(ncname,'a2x6h_Faxa_rainl','standard_name','large_scale_precipitation_flux')
        ncwriteatt(ncname,'a2x6h_Faxa_rainl','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Faxa_rainl','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Faxa_snowc','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_snowc',snowc_assign)
        ncwriteatt(ncname,'a2x6h_Faxa_snowc','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Faxa_snowc','units','kg m-2 s-1')
        ncwriteatt(ncname,'a2x6h_Faxa_snowc','long_name','Convective snow rate (water equivalent)')
        ncwriteatt(ncname,'a2x6h_Faxa_snowc','standard_name','convective_snowfall_flux')
        ncwriteatt(ncname,'a2x6h_Faxa_snowc','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Faxa_snowc','cell_methods','time: mean')

        nccreate(ncname,'a2x6h_Faxa_snowl','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_snowl',snowl_assign)
        ncwriteatt(ncname,'a2x6h_Faxa_snowl','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Faxa_snowl','units','kg m-2 s-1')
        ncwriteatt(ncname,'a2x6h_Faxa_snowl','long_name','Large-scale (stable) snow rate (water equivalent)')
        ncwriteatt(ncname,'a2x6h_Faxa_snowl','standard_name','large_scale_snowfall_flux')
        ncwriteatt(ncname,'a2x6h_Faxa_snowl','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Faxa_snowl','cell_methods','time: mean')


        nccreate(ncname,'a2x6h_Sa_u','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_u',u_assign)
        ncwriteatt(ncname,'a2x6h_Sa_u','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Sa_u','units','m s-1')
        ncwriteatt(ncname,'a2x6h_Sa_u','long_name','Zonal wind at the lowest model level')
        ncwriteatt(ncname,'a2x6h_Sa_u','standard_name','eastward_wind')
        ncwriteatt(ncname,'a2x6h_Sa_u','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Sa_u','cell_methods','time: mean')


        nccreate(ncname,'a2x6h_Sa_v','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_v',v_assign)
        ncwriteatt(ncname,'a2x6h_Sa_v','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x6h_Sa_v','units','m s-1')
        ncwriteatt(ncname,'a2x6h_Sa_v','long_name','Meridional wind at the lowest model level')
        ncwriteatt(ncname,'a2x6h_Sa_v','standard_name','northward_wind')
        ncwriteatt(ncname,'a2x6h_Sa_v','internal_dname','a2x6h')
        ncwriteatt(ncname,'a2x6h_Sa_v','cell_methods','time: mean')

        ncwriteatt(ncname,'/','creation_method', 'Scaling CAM4 regional product (ds199.1) to tower site level met')
        ncwriteatt(ncname,'/','creation_date', datestr(now))
        ncwriteatt(ncname,'/','author', 'Brett Raczka, bmraczka@ucar.edu')
        ncwriteatt(ncname,'/','code', '/glade/work/bmraczka/SIFMIP2/CAM4_biascorrect_SIFMIP2.m')
    
       end %% ensemble assignment loop  

% Important to clear the CAM4 file variables because leap year makes file size inconsistent

clear Faxa_swndr Faxa_swvdr Faxa_swndf Faxa_swvdf Faxa_rainc Faxa_rainl Faxa_snowc Faxa_snowl 
clear  Sa_u Sa_v Sa_tbot Sa_shum Sa_pbot Faxa_lwdn Sa_z Sa_ptem Sa_dens Sa_pslv  

end   %% Year Loop

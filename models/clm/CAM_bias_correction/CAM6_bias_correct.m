clear all
close all

% This script applies a bias correction to CAM6 reanalysis meteorology roughly
% following Yong-Fei Zhang (UT-Austin) dissertation approach. Meteorology data
% from the site location of interest should be used as 'truth'.
% All met variables use 'scaling' approach with exception for snow/rain precip.
% This script is intended to be used after running CAM6_site_grid.sh 

% Enter in site location data. US-NR1 flux tower used as example
 SITE_lat=40.03;
 SITE_lon=-105.55+360; % degrees East
 SITE_doma_area=0.000274675326661479; % radian^2 cell area (unique to grid location)
 SITE_doma_aream=0.000274672687732932; % radian^2 cell area mapping (same value as above)
 SITE_doma_mask=1   ; % all values =1 
 SITE_doma_frac=1   ; % all values =1 area fraction


% Input data files
 path_CAM_forcing = '/glade/work/bmraczka/CAM6_NR1/';
 path_towermet_forcing = '/glade/work/bmraczka/SIFMIP2/tower_met_forcing/';
% Output data files
 path_scaled_CAM_forcing = '/glade/work/bmraczka/CAM6_NR1/scaled/';


% Include Diagnostics?  1=ON, 2=OFF
Diagnostics=1;


% Site level grid extract CAM6 naming convention:  

% CAM6_NR1.cpl_NINST.ha2x1hi.YEAR.nc   % SOLAR:     a2x1hi_Faxa_swndr,a2x1hi_Faxa_swvdr,a2x1hi_Faxa_swndf,a2x1hi_Faxa_swvdf   
% CAM6_NR1.cpl_NINST.ha2x3h.YEAR.nc    % NON-SOLAR: a2x3h_Faxa_rainc,a2x3h_Faxa_rainl,a2x3h_Faxa_snowc,a2x3h_Faxa_snowl,a2x3h_Faxa_lwdn
% CAM6_NR1.cpl_NINST.ha2x1h.YEAR.nc    % 1hr state: a2x1h_Sa_u,a2x1h_Sa_v 
% CAM6_NR1.cpl_NINST.ha2x3h.YEAR.nc    % 3hr state; a2x3h_Sa_tbot,a2x3h_Sa_shum,a2x3h_Sa_pbot
% CAM6_NR1.cpl_NINST.ha2x3h.YEAR.nc    % unchanged 3hr state: Sa_z Sa_ptem Sa_dens Sa_pslv Sa_topo


% 'yearstr' is list of years where both CAM6 reanalysis and tower met forcing
% is available and the years we will create the bias-correct/scaled reanalysis product
yearstr={'2011','2012','2013','2014','2015','2016','2017','2018','2019'};
% 'yeartower' adds previous year to cover the UTC to MST shift 
yeartower={'2010','2011','2012','2013','2014','2015','2016','2017','2018','2019'};


enstr={'0001','0002','0003','0004','0005','0006','0007','0008','0009','0010',...
         '0011','0012','0013','0014','0015','0016','0017','0018','0019','0020',...
         '0021','0022','0023','0024','0025','0026','0027','0028','0029','0030',...
         '0031','0032','0033','0034','0035','0036','0037','0038','0039','0040',...
         '0041','0042','0043','0044','0045','0046','0047','0048','0049','0050',...
         '0051','0052','0053','0054','0055','0056','0057','0058','0059','0060',...
         '0061','0062','0063','0064','0065','0066','0067','0068','0069','0070',...
         '0071','0072','0073','0074','0075','0076','0077','0078','0079','0080'};



% Site level met forcing uses PLUMBER2 protocol 
% time -->30 min increments in LST (MST)
% Tair-->Kelvin
% Qair--> specific humidity (kg kg-1)
% Wind--> (m/s)
% SWdown --> (W/m2)
% LWdown ---> (W/m2)
% Precip --> (mm/s) or (kg m-2 s-1)
% Psurf  --> (Pa)


% Uploading all years of tower forcing, then selecting a particular year
TBOT_master=[];
SH_master=[];
WIND_master=[];
FSDS_master=[];
FLDS_master=[];
PRECTmms_master=[];
PSRF_master=[];
YEAR_master=[];

for ii = 1:length(yeartower);
 
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

% Year Loop
for ii = 1:length(yearstr);

% Time Zones:
% Tower Meterology Forcing (PLUMBER):  LST (MST):  UTC-7
% CAM6 Reanalsysis: UTC

% Syncing time steps:
% Tower Forcing (LST) must be advanced 7 hours to synchronize with CAM6 (UTC)
% In addition there is a 6 hour forward shift in the CAM6 reanalysis yearly files
% Total shift of Tower Forcing: +13 hours with exception of
% SW radiation which was shifted +12 hours for better synchronization

% Both CAM6 and Tower forcing include leap days


% Select only indices of YEAR ii
indices=find(YEAR_master==str2num(yearstr{ii}));

% Pushing met forcing 13 hours forward (go back 26 indices)
ta_30=TBOT_master(indices(1)-26:indices(end)-26)';                        % TBOT (Kelvin)
q_30=SH_master(indices(1)-26:indices(end)-26)';                           % SH (kg kg-1)
wind_30=WIND_master(indices(1)-26:indices(end)-26)';                      % Total Wind (m/s) 
sw_30=FSDS_master(indices(1)-24:indices(end)-24)';                        % FSDS (W/m2)
lw_30=FLDS_master(indices(1)-26:indices(end)-26)';                        % FLDS (W/m2)
ppt_30=PRECTmms_master(indices(1)-26:indices(end)-26)';                   % PRECTmms  (mm/s) or (kg m-2 s-1)
ps_30=PSRF_master(indices(1)-26:indices(end)-26)';                        % PSRF (Pa)

clear indices

% Generate 1 and 3 hourly averages from tower met forcing to compare against
% CAM6 reanalysis 

ta_1hr=mean(reshape(ta_30,2,numel(ta_30)/2),1);        
q_1hr=mean(reshape(q_30,2,numel(q_30)/2),1);           
wind_1hr=mean(reshape(wind_30,2,numel(wind_30)/2),1);  
sw_1hr=mean(reshape(sw_30,2,numel(sw_30)/2),1);        
lw_1hr=mean(reshape(lw_30,2,numel(lw_30)/2),1);        
ppt_1hr=mean(reshape(ppt_30,2,numel(ppt_30)/2),1);     
ps_1hr=mean(reshape(ps_30,2,numel(ps_30)/2),1);


ta_3hr=mean(reshape(ta_1hr,3,numel(ta_1hr)/3),1);       
q_3hr=mean(reshape(q_1hr,3,numel(q_1hr)/3),1);         
wind_3hr=mean(reshape(wind_1hr,3,numel(wind_1hr)/3),1); 
sw_3hr=mean(reshape(sw_1hr,3,numel(sw_1hr)/3),1);       
lw_3hr=mean(reshape(lw_1hr,3,numel(lw_1hr)/3),1);       
ppt_3hr=mean(reshape(ppt_1hr,3,numel(ppt_1hr)/3),1);    
ps_3hr=mean(reshape(ps_1hr,3,numel(ps_1hr)/3),1);       

clear ta_30 q_30 wind_30 sw_30 lw_30 ppt_30 ps_30


   % Ensemble Loop 
   for jj = 1:80;  

       %Downloading all ensemble members for all  reanalysis variables
    
       %SOLAR
       Faxa_swndr(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x1hi.' yearstr{ii} '.nc'],'a2x1hi_Faxa_swndr'); % 1 hour resolution
       Faxa_swvdr(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x1hi.' yearstr{ii} '.nc'],'a2x1hi_Faxa_swvdr'); % 1 hour resolution
       Faxa_swndf(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x1hi.' yearstr{ii} '.nc'],'a2x1hi_Faxa_swndf'); % 1 hour resolution
       Faxa_swvdf(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x1hi.' yearstr{ii} '.nc'],'a2x1hi_Faxa_swvdf'); % 1 hour resolution
    
       %NONSOLAR
       Faxa_rainc(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Faxa_rainc'); % 3 hour resolution
       Faxa_rainl(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Faxa_rainl'); % 3 hour resolution
       Faxa_snowc(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Faxa_snowc'); % 3 hour resolution
       Faxa_snowl(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Faxa_snowl'); % 3 hour resolution

       % 1hr STATE
       Sa_u(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x1h.' yearstr{ii} '.nc'],'a2x1h_Sa_u'); % 1 hour resolution
       Sa_v(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x1h.' yearstr{ii} '.nc'],'a2x1h_Sa_v'); % 1 hour resolution
    
       % 3hr STATE
       Sa_tbot(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Sa_tbot'); % 3 hour resolution
       Sa_shum(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Sa_shum'); % 3 hour resolution
       Sa_pbot(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Sa_pbot');       % 3 hour resolution
       Faxa_lwdn(jj,:)  = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Faxa_lwdn');  % 3 hour resolution
       
       % (3hr STATE Unchanged)
       Sa_z(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Sa_z');       % 3 hour resolution
       Sa_ptem(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Sa_ptem'); % 3 hour resolution
       Sa_dens(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Sa_dens'); % 3 hour resolution
       Sa_pslv(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Sa_pslv'); % 3 hour resolution
       Sa_topo(jj,:) = ncread([path_CAM_forcing enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'],'a2x3h_Sa_topo'); % 3 hour resolution
       

   end %% ensemble loop
     
      %% Initialize the ensemble mean scale values to the default of '0' which means no change
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
       Sa_topo_scale =zeros(1,length(Sa_topo(1,:)));
      
       % Calculate the ensemble mean for each CAM6 variable
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

        % Calculate CAM6 total wind vector for comparison with tower forcing
        Sa_wind_mean= (Sa_v_mean.^2+Sa_u_mean.^2).^0.5;

       Sa_tbot_mean = mean(Sa_tbot,1);
       Sa_shum_mean = mean(Sa_shum,1);
       Sa_pbot_mean = mean(Sa_pbot,1);
       Faxa_lwdn_mean = mean(Faxa_lwdn,1);


     % Need this initialization for precip loop
       Faxa_rain=Faxa_rainl;
       Faxa_snow=Faxa_snowl;


      %% Calculate the scaled correction factors for each year 

      %% SOLAR 1 hour resolution
      for ind=1:(length(sw_1hr));

        if (Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind)+Faxa_swndf_mean(ind)+Faxa_swvdr_mean(ind)>0)      
      
           Faxa_swndr_scale(ind)= Faxa_swndr_mean(ind) .* (sw_1hr(ind)./(Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind) ...
                                                             +Faxa_swndf_mean(ind)+Faxa_swvdf_mean(ind)));
           Faxa_swvdr_scale(ind)= Faxa_swvdr_mean(ind) .* (sw_1hr(ind)./(Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind) ...
                                                             +Faxa_swndf_mean(ind)+Faxa_swvdf_mean(ind)));
           Faxa_swndf_scale(ind)= Faxa_swndf_mean(ind) .* (sw_1hr(ind)./(Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind) ...
                                                             +Faxa_swndf_mean(ind)+Faxa_swvdf_mean(ind)));
           Faxa_swvdf_scale(ind)= Faxa_swvdf_mean(ind) .* (sw_1hr(ind)./(Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind) ...
                                                             +Faxa_swndf_mean(ind)+Faxa_swvdf_mean(ind)));
        end
      end 

      %% NON-SOLAR 3 hour resolution      
      for ind=1:(length(ppt_3hr));
      
      % Precip loop needs to query adjusted temperature to re-assign precip to snow or rain

          if ta_3hr(ind)>273.15
             % Forcing all unadjusted precip ensemble members to new variable Faxa_rain according to adjusted temperature
             Faxa_rain(:,ind)=Faxa_rainc(:,ind)+Faxa_rainl(:,ind)+Faxa_snowc(:,ind)+Faxa_snowl(:,ind);
             Faxa_snow(:,ind)=0;

             % Adjustment will only be applied to total rain. CLM does not care if precip is convective/large scale 
             Faxa_rainc_scale(ind)= 0;
             Faxa_rainl_scale(ind)= 0;
             Faxa_snowc_scale(ind)= 0;
             Faxa_snowl_scale(ind)= 0;
             Faxa_rain_scale(ind)= ppt_3hr(ind);
          else
            % Forcing all unadjusted precip ensemble members to new variable Faxa_ snow according to adjusted temperature
            Faxa_snow(:,ind)= Faxa_rainc(:,ind)+Faxa_rainl(:,ind)+Faxa_snowc(:,ind)+Faxa_snowl(:,ind);
            Faxa_rain(:,ind)=0;

            % Adjustment will only be applied to total snow. CLM does not care if precip is convective/large scale
             Faxa_rainc_scale(ind)=0;
             Faxa_rainl_scale(ind)=0;
             Faxa_snowc_scale(ind)=0;
             Faxa_snowl_scale(ind)=0;
             Faxa_snow_scale(ind)= ppt_3hr(ind);
          end

      end

       Faxa_rain_mean= mean(Faxa_rain,1);
       Faxa_snow_mean= mean(Faxa_snow,1);

       % 1hr STATE
       Sa_u_scale = Sa_u_mean.*(wind_1hr./Sa_wind_mean); 
       Sa_v_scale = Sa_v_mean.*(wind_1hr./Sa_wind_mean);

       % 3hr STATE
       Sa_tbot_scale = Sa_tbot_mean.*(ta_3hr./Sa_tbot_mean);
       Sa_shum_scale = Sa_shum_mean.*(q_3hr./Sa_shum_mean);
       Sa_pbot_scale = Sa_pbot_mean.*(ps_3hr./Sa_pbot_mean);
       Faxa_lwdn_scale = Faxa_lwdn_mean.*(lw_3hr./Faxa_lwdn_mean);

       % 3hr STATE (Unchanged - no scaling)
       %Sa_z_scale    =  Sa_z_mean;
       %Sa_ptem_scale =  Sa_ptem_mean;
       %Sa_dens_scale =  Sa_desn_mean;
       %Sa_pslv_scale =  Sa_pslv_mean;
       %Sa_topo_scale =  Sa_topo_mean;

   %% Apply the scaled corrections identically to each ensemble member, then generate corrected CAM6 reanalysis
   %% ensemble_scaled=ensemble + (ensemble_mean_scaled-ensemble_mean)
   %% Put clamping (hard bounds) on non-physical values.

       % SOLAR
       % e.g Faxa_swndr(ensemble(80),time(1hr,3hr)
       Faxa_swndr_adjust  = Faxa_swndr+repmat(Faxa_swndr_scale-Faxa_swndr_mean,80,1); Faxa_swndr_adjust(Faxa_swndr_adjust<0)=0;  
       Faxa_swvdr_adjust  = Faxa_swvdr+repmat(Faxa_swvdr_scale-Faxa_swvdr_mean,80,1); Faxa_swvdr_adjust(Faxa_swvdr_adjust<0)=0;
       Faxa_swndf_adjust  = Faxa_swndf+repmat(Faxa_swndf_scale-Faxa_swndf_mean,80,1); Faxa_swndf_adjust(Faxa_swndf_adjust<0)=0;
       Faxa_swvdf_adjust  = Faxa_swvdf+repmat(Faxa_swvdf_scale-Faxa_swvdf_mean,80,1); Faxa_swvdf_adjust(Faxa_swvdf_adjust<0)=0;

       %NONSOLAR
       % Purposely setting convective rain and snow to zero
       % Purposely setting snow and rain adjusted variables to large scale snow (snowl) and rain (rainl)
       % These will carry the ensemble spread for snow/rain
 
       Faxa_rainc_adjust = Faxa_rainc+repmat(Faxa_rainc_scale-Faxa_rainc_mean,80,1); Faxa_rainc_adjust(:,:)=0;
       Faxa_rainl_adjust = Faxa_rain+repmat(Faxa_rain_scale-Faxa_rain_mean,80,1); Faxa_rainl_adjust(Faxa_rainl_adjust<0)=0;
       Faxa_snowc_adjust = Faxa_snowc+repmat(Faxa_snowc_scale-Faxa_snowc_mean,80,1); Faxa_snowc_adjust(:,:)=0;
       Faxa_snowl_adjust = Faxa_snow+repmat(Faxa_snow_scale-Faxa_snow_mean,80,1); Faxa_snowl_adjust(Faxa_snowl_adjust<0)=0;

       % 1hr STATE, meridional and zonal winds allowed to go negative
       Sa_u_adjust = Sa_u+repmat(Sa_u_scale-Sa_u_mean,80,1);
       Sa_v_adjust = Sa_v+repmat(Sa_v_scale-Sa_v_mean,80,1); 

       % 3hr STATE
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
if Diagnostics==1

     % Diagnostic Figures Only
     % Comparing original CAM6 and Tower met forcing
     % Check for time synchronization, ensemble adjustment, snow/rain partitioning etc. 
     
       figure(1)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       j=plot(mean(Faxa_swndr,1)+ ...
              mean(Faxa_swndf,1)+ ...
              mean(Faxa_swvdr,1)+ ...
              mean(Faxa_swvdf,1), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(sw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Total Shortwave (W m^{-2})','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr{ii}] , 'FontSize',fontpt);

       axes(ha(2));
       j=plot(repelem(mean(Faxa_rainl,1),1,3)'+ ...
              repelem(mean(Faxa_rainc,1),1,3)'+ ...
              repelem(mean(Faxa_snowl,1),1,3)'+ ...
              repelem(mean(Faxa_snowc,1),1,3)', ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ppt_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)

       axes(ha(3));
       j=plot(repelem(mean(Sa_tbot,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ta_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Temperature (K)','FontSize', fontpt)

       axes(ha(4));
       j=plot(repelem(mean(Sa_shum,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(q_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

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

       j=plot(mean(Sa_u,1), ...
              '--', 'color',rgb('darksalmon'), 'linewidth', mean_width); hold on
       k=plot(mean(Sa_v,1), ...
              '-', 'color',rgb('darksalmon'), 'linewidth', mean_width);
       l=plot( (mean(Sa_u,1).^2+ ...
               mean(Sa_v,1).^2).^0.5, ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width);

       m=plot(wind_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k,l,m], 'Zonal, CAM6 Ens Mean', 'Meridional, CAM6 Ens Mean', ...
                     'Total, CAM6 Ens Mean','Tower Met','FontSize',fontpt);

       ylabel('Wind (m s^{-1})','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr{ii}] , 'FontSize',fontpt);


       axes(ha(2));
       j=plot(repelem(mean(Sa_tbot,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ta_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Temperature (K)','FontSize', fontpt)

       axes(ha(3));
       j=plot(repelem(mean(Faxa_rainl,1),1,3)'+ ...
              repelem(mean(Faxa_rainc,1),1,3)'+ ...
              repelem(mean(Faxa_snowl,1),1,3)'+ ...
              repelem(mean(Faxa_snowc,1),1,3)', ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ppt_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

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
       j=plot(repelem(mean(Sa_pbot,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ps_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('Pressure (Pa)','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr{ii}] , 'FontSize',fontpt);

       axes(ha(2));
       j=plot(repelem(mean(Faxa_lwdn,1),1,3), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(lw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);

       ylabel('LW radiation (W m^{-1})','FontSize', fontpt)
       linkaxes(ha(1:2), 'x');

display('  ')
display('Finished tower CAM6 diagnostic plots. Press enter to proceed to adjusted CAM6 plots..')
display('  ')
pause
close all


     % For plotting purposes only, make these 3 hour products into 1 hour products for easier comparison 
     % Matrices (80, 2920) --> (80, 8760); 

       figure(4)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3)/2 scrsz(4)]);
       fontpt=16;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       ax1=plot(Faxa_swndr(:,188*24:195*24)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndr_adjust(:,188*24:195*24)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width); 
       j=plot(mean(Faxa_swndr(:,188*24:195*24),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); 
       k=plot(mean(Faxa_swndr_adjust(:,188*24:195*24),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on 
       legend( [j,k], 'CAM6 Mean', 'Adjusted CAM6 Mean','FontSize',fontpt);       
       ylabel('Near-IR, Direct (W  m^{-2})','FontSize', fontpt)
       title(['CAM6 scaling against tower NR1 Met: July ', yearstr{ii}] , 'FontSize',fontpt);
       xticks([24:24:24*18]);
       set(ha(1),'xticklabel',[]) 

       axes(ha(2));

       ax2=plot(repelem(Faxa_rain(:,188*8:195*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Faxa_rainl_adjust(:,188*8:195*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width); 
       plot(repelem(mean(Faxa_rain(:,188*8:195*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Faxa_rainl_adjust(:,188*8:195*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Rain (mm  s^{-1})','FontSize', fontpt)
       xticks([24:24:24*18]);
       set(ha(2),'xticklabel',[])

       axes(ha(3));

       ax3=plot(repelem(Sa_tbot(:,188*8:195*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_tbot_adjust(:,188*8:195*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width); 
       plot(repelem(mean(Sa_tbot(:,188*8:195*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); 
       plot(repelem(mean(Sa_tbot_adjust(:,188*8:195*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Temperature (K)','FontSize', fontpt)
       xticks([24:24:24*18]);
       set(ha(3),'xticklabel',[])

       axes(ha(4));

       ax4=plot(repelem(Sa_shum(:,188*8:195*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_shum_adjust(:,188*8:195*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Sa_shum(:,188*8:195*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Sa_shum_adjust(:,188*8:195*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Specific Humidity (kg  kg^{-2})','FontSize', fontpt)
       xticks([24:24:24*18]);
       linkaxes(ha(1:4), 'x');

     
 
       figure(5)
       scrsz = get(0,'ScreenSize');
       set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);
       fontpt=12;
       mean_width=1.5;
       member_width=0.5;
       window=10;

       ha=tight_subplot(4,1, [.02,.08],[.10,.06],[.08,.06]);
       axes(ha(1));

       plot(Faxa_swndr(:,32*24:45*24)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndr_adjust(:,32*24:45*24)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       j=plot(mean(Faxa_swndr(:,32*24:45*24),1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); 
       k=plot(mean(Faxa_swndr_adjust(:,32*24:45*24),1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM6 Mean', 'Adjusted CAM6 Mean','FontSize',fontpt); 
       ylabel('Near-IR, Direct (W  m^{-2})','FontSize', fontpt)
       title(['CAM6 scaling against tower NR1 Met: February ', yearstr{ii}] , 'FontSize',fontpt);
       xticks([24:24:24*13]);       

       axes(ha(2));

       plot(repelem(Faxa_snow(:,32*8:45*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Faxa_snowl_adjust(:,32*8:45*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width); 
       plot(repelem(mean(Faxa_snow(:,32*8:45*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); 
       plot(repelem(mean(Faxa_snowl_adjust(:,32*8:45*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Snow (mm  s^{-1})','FontSize', fontpt)
       xticks([24:24:24*13]);

       axes(ha(3));

       plot(repelem(Sa_tbot(:,32*8:45*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_tbot_adjust(:,32*8:45*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width); 
       plot(repelem(mean(Sa_tbot(:,32*8:45*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width); 
       plot(repelem(mean(Sa_tbot_adjust(:,32*8:45*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Temperature (K)','FontSize', fontpt)
       xticks([24:24:24*13]);
       axes(ha(4));

       plot(repelem(Sa_shum(:,32*8:45*8),1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_shum_adjust(:,32*8:45*8),1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Sa_shum(:,32*8:45*8),1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Sa_shum_adjust(:,32*8:45*8),1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Specific Humidity (kg  kg^{-2})','FontSize', fontpt)
       xticks([24:24:24*13]);

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
       k=plot(sw_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);
       title(['CAM6 scaling against tower NR1 Met ', yearstr{ii}] , 'FontSize',fontpt);
       ylabel('Total Shortwave (W m^{-2})','FontSize', fontpt)


       axes(ha(2));

       ax1=plot(Faxa_swndr', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Faxa_swndr_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       j=plot(mean(Faxa_swndr,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       k=plot(mean(Faxa_swndr_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM6 Mean', 'Adjusted CAM6 Mean','FontSize',fontpt);
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
       k=plot(ppt_3hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);
       title(['CAM6 scaling against tower NR1 Met ', yearstr{ii}] , 'FontSize',fontpt);
       ylabel('Total Precipitation (mm s^{-1})','FontSize', fontpt)
       ylim([0 2*10^-3])

       axes(ha(2));
       j=plot(mean(Sa_tbot,1), ...
              '-', 'color',rgb('darkred'), 'linewidth', mean_width); hold on
       k=plot(ta_3hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k], 'CAM6 Ensemble Mean', 'Tower Met','FontSize',fontpt);
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

       m=plot(wind_1hr, '-', 'color',rgb('black'),'linewidth', mean_width);
       grid on
       legend([j,k,l,m], 'Zonal, CAM6 Ens Mean', 'Meridional, CAM6 Ens Mean', ...
                     'Total, CAM6 Ens Mean','Tower Met','FontSize',fontpt);

       ylabel('Wind (m s^{-1})','FontSize', fontpt)
       title(['Niwot Ridge Meteorology ', yearstr{ii}] , 'FontSize',fontpt);
       axes(ha(2));

       ax1=plot(Sa_u', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_u_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);

       j=plot(mean(Sa_u,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       k=plot(mean(Sa_u_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       legend( [j,k], 'CAM6 Mean', 'Adjusted CAM6 Mean','FontSize',fontpt);
       ylabel('Zonal wind (m s^{-1})','FontSize', fontpt)
       axes(ha(3));

       ax2=plot(Sa_v', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(Sa_v_adjust', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(mean(Sa_v,1)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(mean(Sa_v_adjust,1)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Meridional wind (m s^{-1})','FontSize', fontpt)
       axes(ha(4));

       ax3=plot(repelem(Sa_pbot,1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Sa_pbot_adjust,1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Sa_pbot,1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Sa_pbot_adjust,1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Pressure (Pa)','FontSize', fontpt)

       axes(ha(5));

       ax4=plot(repelem(Faxa_lwdn,1,3)', '-', 'color',rgb('lightsalmon'), 'linewidth', member_width); hold on
       plot(repelem(Faxa_lwdn_adjust,1,3)', '-', 'color',rgb('lightgreen'), 'linewidth', member_width);
       plot(repelem(mean(Faxa_lwdn,1),1,3)', '-', 'color',rgb('darkred'), 'linewidth', mean_width);
       plot(repelem(mean(Faxa_lwdn_adjust,1),1,3)', '-', 'color',rgb('forestgreen'), 'linewidth', mean_width);
       grid on
       ylabel('Longwave Radiation (W m^{-2})','FontSize', fontpt)
       linkaxes(ha(1:5), 'x');

display('   ')
display(['Adjusted CAM6 plots and diagnostics completed.  Press enter to write ' yearstr{ii} ' CAM6 adjusted files to netcdf'])
display('   ')
pause
close all


end % Diagnostics Loop
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      % Assign the adjusted CAM variables to netcdf files for each ensemble member/ year
          % CAM6_NR1.cpl_NINST.ha2x1hi.YEAR.nc   % SOLAR:     a2x1hi_Faxa_swndr,a2x1hi_Faxa_swvdr,a2x1hi_Faxa_swndf,a2x1hi_Faxa_swvdf   
          % CAM6_NR1.cpl_NINST.ha2x3h.YEAR.nc    % NON-SOLAR: a2x3h_Faxa_rainc,a2x3h_Faxa_rainl,a2x3h_Faxa_snowc,a2x3h_Faxa_snowl,a2x3h_Faxa_lwdn
          % CAM6_NR1.cpl_NINST.ha2x1h.YEAR.nc    % 1hr state: a2x1h_Sa_u,a2x1h_Sa_v 
          % CAM6_NR1.cpl_NINST.ha2x3h.YEAR.nc    % 3hr state; a2x3h_Sa_tbot,a2x3h_Sa_shum,a2x3h_Sa_pbot
          % CAM6_NR1.cpl_NINST.ha2x3h.YEAR.nc    % unchanged 3hr state: Sa_z Sa_ptem Sa_dens Sa_pslv Sa_topo


          time_SOLAR = ncread([path_CAM_forcing enstr{1} '/CAM6_NR1.cpl_' enstr{1} '.ha2x1hi.' yearstr{ii} '.nc'],'time');
          time_NONSOLAR = ncread([path_CAM_forcing enstr{1} '/CAM6_NR1.cpl_' enstr{1} '.ha2x3h.' yearstr{ii} '.nc'],'time');
          time_1hr = ncread([path_CAM_forcing enstr{1} '/CAM6_NR1.cpl_' enstr{1} '.ha2x1h.' yearstr{ii} '.nc'],'time');
          time_3hr= ncread([path_CAM_forcing enstr{1} '/CAM6_NR1.cpl_' enstr{1} '.ha2x3h.' yearstr{ii} '.nc'],'time');

          dim_time_SOLAR=length(time_SOLAR);        % ha2x1hi
          dim_time_NONSOLAR=length(time_NONSOLAR);  % ha2x3h
          dim_time_1hr=length(time_1hr);            % ha2x1h
          dim_time_3hr=length(time_3hr);            % ha2x3h

          time_bnds_SOLAR = ncread([path_CAM_forcing enstr{1} '/CAM6_NR1.cpl_' enstr{1} '.ha2x1hi.' yearstr{ii} '.nc'],'time_bnds');
          time_bnds_NONSOLAR = ncread([path_CAM_forcing enstr{1} '/CAM6_NR1.cpl_' enstr{1} '.ha2x3h.' yearstr{ii} '.nc'],'time_bnds');
          time_bnds_1hr = ncread([path_CAM_forcing enstr{1} '/CAM6_NR1.cpl_' enstr{1} '.ha2x1h.' yearstr{ii} '.nc'],'time_bnds');
          time_bnds_3hr = ncread([path_CAM_forcing enstr{1} '/CAM6_NR1.cpl_' enstr{1} '.ha2x3h.' yearstr{ii} '.nc'],'time_bnds');

       
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
        topo_assign=ones(1,1,length(Sa_topo(1,:)))*NaN;

        for jj = 1:80;  %% ensemble loop         
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
        topo_assign(1,1,:)=Sa_topo(jj,:);

        % Write SOLAR (ha2x1hi) File

        ncname=[path_scaled_CAM_forcing  enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x1hi.' yearstr{ii} '.nc'];
        ncid=netcdf.create(ncname,'CLOBBER');
        netcdf.close(ncid);
        nccreate(ncname,'a2x1hi_Faxa_swndf','Dimensions',{'a2x1hi_nx',1,'a2x1hi_ny',1,'time',Inf},'Datatype','single','Format','classic'); 
        ncwrite(ncname,'a2x1hi_Faxa_swndf',swndf_assign)
        ncwriteatt(ncname,'a2x1hi_Faxa_swndf','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x1hi_Faxa_swndf','units','W m-2')
        ncwriteatt(ncname,'a2x1hi_Faxa_swndf','long_name','Diffuse near-infrared incident solar radiation')
        ncwriteatt(ncname,'a2x1hi_Faxa_swndf','standard_name','surface_downward_diffuse_shortwave_flux_due_to_near_infrared_radiation')
        ncwriteatt(ncname,'a2x1hi_Faxa_swndf','internal_dname','a2x1hi')

        nccreate(ncname,'a2x1hi_Faxa_swndr','Dimensions',{'a2x1hi_nx',1,'a2x1hi_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x1hi_Faxa_swndr',swndr_assign)
        ncwriteatt(ncname,'a2x1hi_Faxa_swndr','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x1hi_Faxa_swndr','units','W m-2')
        ncwriteatt(ncname,'a2x1hi_Faxa_swndr','long_name','Direct near-infrared incident solar radiation')
        ncwriteatt(ncname,'a2x1hi_Faxa_swndr','standard_name','surface_downward_direct_shortwave_flux_due_to_near_infrared_radiation')
        ncwriteatt(ncname,'a2x1hi_Faxa_swndr','internal_dname','a2x1hi')
               
        nccreate(ncname,'a2x1hi_Faxa_swvdf','Dimensions',{'a2x1hi_nx',1,'a2x1hi_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x1hi_Faxa_swvdf',swvdf_assign)
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdf','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdf','units','W m-2')
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdf','long_name','Diffuse visible incident solar radiation')
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdf','standard_name','surface_downward_diffuse_shortwave_flux_due_to_visible_radiation')
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdf','internal_dname','a2x1hi')
    
        nccreate(ncname,'a2x1hi_Faxa_swvdr','Dimensions',{'a2x1hi_nx',1,'a2x1hi_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x1hi_Faxa_swvdr',swvdr_assign)
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdr','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdr','units','W m-2')
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdr','long_name','Direct visible incident solar radiation')
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdr','standard_name','surface_downward_direct_shortwave_flux_due_to_visible_radiation')
        ncwriteatt(ncname,'a2x1hi_Faxa_swvdr','internal_dname','a2x1hi')

        nccreate(ncname,'doma_area','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_area',SITE_doma_area)
        ncwriteatt(ncname,'doma_area','_FillValue', 0.)
        ncwriteatt(ncname,'doma_area','units','radian^2')
        ncwriteatt(ncname,'doma_area','long_name','cell_area_model')
        ncwriteatt(ncname,'doma_area','standard_name','cell area from model')
        ncwriteatt(ncname,'doma_area','internal_dname','doma')

        nccreate(ncname,'doma_aream','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_aream',SITE_doma_aream)
        ncwriteatt(ncname,'doma_aream','_FillValue', 0.)
        ncwriteatt(ncname,'doma_aream','units','radian^2')
        ncwriteatt(ncname,'doma_aream','long_name','cell_area_mapping')
        ncwriteatt(ncname,'doma_aream','standard_name','cell area from mapping file')
        ncwriteatt(ncname,'doma_aream','internal_dname','doma')

        nccreate(ncname,'doma_frac','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_frac',SITE_doma_frac)
        ncwriteatt(ncname,'doma_frac','_FillValue', 0.) 
        ncwriteatt(ncname,'doma_frac','units','1')
        ncwriteatt(ncname,'doma_frac','long_name','area_fraction')
        ncwriteatt(ncname,'doma_frac','standard_name','area fraction')
        ncwriteatt(ncname,'doma_frac','internal_dname','doma')

        nccreate(ncname,'doma_lat','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lat',SITE_lat)
        ncwriteatt(ncname,'doma_lat','_FillValue', 0.)
        ncwriteatt(ncname,'doma_lat','units','degrees north')
        ncwriteatt(ncname,'doma_lat','long_name','latitude')
        ncwriteatt(ncname,'doma_lat','standard_name','latitude')
        ncwriteatt(ncname,'doma_lat','internal_dname','doma')
        
        nccreate(ncname,'doma_lon','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lon',SITE_lon)
        ncwriteatt(ncname,'doma_lon','_FillValue', 0.)
        ncwriteatt(ncname,'doma_lon','units','degrees east')
        ncwriteatt(ncname,'doma_lon','long_name','longitude')
        ncwriteatt(ncname,'doma_lon','standard_name','longitude')
        ncwriteatt(ncname,'doma_lon','internal_dname','doma')

        nccreate(ncname,'doma_mask','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_mask',SITE_doma_mask)
        ncwriteatt(ncname,'doma_mask','_FillValue', 0.)
        ncwriteatt(ncname,'doma_mask','units','1')
        ncwriteatt(ncname,'doma_mask','long_name','mask')
        ncwriteatt(ncname,'doma_mask','standard_name','mask')
        ncwriteatt(ncname,'doma_mask','internal_dname','doma')

        nccreate(ncname,'time','Dimensions',{'time',Inf},'Format','classic');
        ncwrite(ncname,'time',time_SOLAR)

        if str2num(yearstr{ii})==2011
        ncwriteatt(ncname,'time','units',['days since 2010-07-01 00:00:00'])
        else
        ncwriteatt(ncname,'time','units',['days since 2011-01-01 00:00:00'])
        end

        ncwriteatt(ncname,'time','calendar','gregorian')
        ncwriteatt(ncname,'time','bounds','time_bnds')

        nccreate(ncname,'time_bnds','Dimensions',{'ntb',2,'time',Inf},'Format','classic');
        ncwrite(ncname,'time_bnds',time_bnds_SOLAR)

        ncwriteatt(ncname,'/','creation_method', 'Scaling CAM6 regional product to tower site level met')
        ncwriteatt(ncname,'/','creation_date', datestr(now))
        ncwriteatt(ncname,'/','author', 'Brett Raczka, bmraczka@ucar.edu')


        % Write NON-SOLAR (ha2x3h) File
        ncname=[path_scaled_CAM_forcing  enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x3h.' yearstr{ii} '.nc'];
        ncid=netcdf.create(ncname,'CLOBBER');
        netcdf.close(ncid);
        nccreate(ncname,'time','Dimensions',{'time',Inf},'Format','classic'); 
        ncwrite(ncname,'time',time_NONSOLAR)

        if str2num(yearstr{ii})==2011
        ncwriteatt(ncname,'time','units',['days since 2010-07-01 00:00:00'])
        else
        ncwriteatt(ncname,'time','units',['days since 2011-01-01 00:00:00'])
        end
        
        ncwriteatt(ncname,'time','calendar','gregorian')
        ncwriteatt(ncname,'time','bounds','time_bnds')
        ncwriteatt(ncname,'time','cell_methods','time: mean')

        nccreate(ncname,'time_bnds','Dimensions',{'ntb',2,'time',Inf},'Format','classic');
        ncwrite(ncname,'time_bnds',time_bnds_NONSOLAR)
        ncwriteatt(ncname,'time_bnds','cell_methods','time: mean')

        nccreate(ncname,'doma_lat','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lat',SITE_lat)
        ncwriteatt(ncname,'doma_lat','_FillValue', 0.)
        ncwriteatt(ncname,'doma_lat','units','degrees north')
        ncwriteatt(ncname,'doma_lat','long_name','latitude')
        ncwriteatt(ncname,'doma_lat','standard_name','latitude')
        ncwriteatt(ncname,'doma_lat','internal_dname','doma')

        nccreate(ncname,'doma_lon','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lon',SITE_lon)
        ncwriteatt(ncname,'doma_lon','_FillValue', 0.)
        ncwriteatt(ncname,'doma_lon','units','degrees east')
        ncwriteatt(ncname,'doma_lon','long_name','longitude')
        ncwriteatt(ncname,'doma_lon','standard_name','longitude')
        ncwriteatt(ncname,'doma_lon','internal_dname','doma')

        nccreate(ncname,'doma_area','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_area',SITE_doma_area)
        ncwriteatt(ncname,'doma_area','_FillValue', 0.)
        ncwriteatt(ncname,'doma_area','units','radian^2')
        ncwriteatt(ncname,'doma_area','long_name','cell_area_model')
        ncwriteatt(ncname,'doma_area','standard_name','cell area from model')
        ncwriteatt(ncname,'doma_area','internal_dname','doma')

        nccreate(ncname,'doma_aream','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_aream',SITE_doma_aream)
        ncwriteatt(ncname,'doma_aream','_FillValue', 0.)
        ncwriteatt(ncname,'doma_aream','units','radian^2')
        ncwriteatt(ncname,'doma_aream','long_name','cell_area_mapping')
        ncwriteatt(ncname,'doma_aream','standard_name','cell area from mapping file')
        ncwriteatt(ncname,'doma_aream','internal_dname','doma')

        nccreate(ncname,'doma_mask','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_mask',SITE_doma_mask)
        ncwriteatt(ncname,'doma_mask','_FillValue', 0.)
        ncwriteatt(ncname,'doma_mask','units','1')
        ncwriteatt(ncname,'doma_mask','long_name','mask')
        ncwriteatt(ncname,'doma_mask','standard_name','mask')
        ncwriteatt(ncname,'doma_mask','internal_dname','doma')

        nccreate(ncname,'doma_frac','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_frac',SITE_doma_frac)
        ncwriteatt(ncname,'doma_frac','_FillValue', 0.)
        ncwriteatt(ncname,'doma_frac','units','1')
        ncwriteatt(ncname,'doma_frac','long_name','area_fraction')
        ncwriteatt(ncname,'doma_frac','standard_name','area fraction')
        ncwriteatt(ncname,'doma_frac','internal_dname','doma')
                                                                                                                        
        nccreate(ncname,'a2x3h_Sa_z','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Sa_z',z_assign)
        ncwriteatt(ncname,'a2x3h_Sa_z','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Sa_z','units','m')
        ncwriteatt(ncname,'a2x3h_Sa_z','long_name','Height at the lowest model level')
        ncwriteatt(ncname,'a2x3h_Sa_z','standard_name','height')
        ncwriteatt(ncname,'a2x3h_Sa_z','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Sa_z','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Sa_topo','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Sa_topo',topo_assign)
        ncwriteatt(ncname,'a2x3h_Sa_topo','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Sa_topo','units','m')
        ncwriteatt(ncname,'a2x3h_Sa_topo','long_name','Surface height')
        ncwriteatt(ncname,'a2x3h_Sa_topo','standard_name','height')
        ncwriteatt(ncname,'a2x3h_Sa_topo','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Sa_topo','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Sa_tbot','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Sa_tbot',tbot_assign)
        ncwriteatt(ncname,'a2x3h_Sa_tbot','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Sa_tbot','units','K')
        ncwriteatt(ncname,'a2x3h_Sa_tbot','long_name','Temperature at the lowest model level')
        ncwriteatt(ncname,'a2x3h_Sa_tbot','standard_name','air_temperature')
        ncwriteatt(ncname,'a2x3h_Sa_tbot','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Sa_tbot','cell_methods','time: mean')
        
        nccreate(ncname,'a2x3h_Sa_ptem','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Sa_ptem',ptem_assign)
        ncwriteatt(ncname,'a2x3h_Sa_ptem','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Sa_ptem','units','K')
        ncwriteatt(ncname,'a2x3h_Sa_ptem','long_name','Potential temperature at the lowest model level')
        ncwriteatt(ncname,'a2x3h_Sa_ptem','standard_name','air_potential_temperature')
        ncwriteatt(ncname,'a2x3h_Sa_ptem','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Sa_ptem','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Sa_shum','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Sa_shum',shum_assign)
        ncwriteatt(ncname,'a2x3h_Sa_shum','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Sa_shum','units','kg kg-1')
        ncwriteatt(ncname,'a2x3h_Sa_shum','long_name','Specific humidity at the lowest model level')
        ncwriteatt(ncname,'a2x3h_Sa_shum','standard_name','specific humidity')
        ncwriteatt(ncname,'a2x3h_Sa_shum','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Sa_shum','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Sa_dens','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Sa_dens',dens_assign)
        ncwriteatt(ncname,'a2x3h_Sa_dens','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Sa_dens','units','kg m-3')
        ncwriteatt(ncname,'a2x3h_Sa_dens','long_name','Density at the lowest model level')
        ncwriteatt(ncname,'a2x3h_Sa_dens','standard_name','air_density')
        ncwriteatt(ncname,'a2x3h_Sa_dens','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Sa_dens','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Sa_pbot','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Sa_pbot',pbot_assign)
        ncwriteatt(ncname,'a2x3h_Sa_pbot','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Sa_pbot','units','Pa')
        ncwriteatt(ncname,'a2x3h_Sa_pbot','long_name','Pressure at the lowest model level')
        ncwriteatt(ncname,'a2x3h_Sa_pbot','standard_name','air_pressure')
        ncwriteatt(ncname,'a2x3h_Sa_pbot','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Sa_pbot','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Sa_pslv','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Sa_pslv',pslv_assign)
        ncwriteatt(ncname,'a2x3h_Sa_pslv','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Sa_pslv','units','Pa')
        ncwriteatt(ncname,'a2x3h_Sa_pslv','long_name','Sea level pressure')
        ncwriteatt(ncname,'a2x3h_Sa_pslv','standard_name','air_pressure_at_sea_level')
        ncwriteatt(ncname,'a2x3h_Sa_pslv','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Sa_pslv','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Faxa_lwdn','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Faxa_lwdn',lwdn_assign)
        ncwriteatt(ncname,'a2x3h_Faxa_lwdn','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Faxa_lwdn','units','W m-2')
        ncwriteatt(ncname,'a2x3h_Faxa_lwdn','long_name','Downward longwave heat flux')
        ncwriteatt(ncname,'a2x3h_Faxa_lwdn','standard_name','downwelling_longwave_flux')
        ncwriteatt(ncname,'a2x3h_Faxa_lwdn','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Faxa_lwdn','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Faxa_rainc','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Faxa_rainc',rainc_assign)
        ncwriteatt(ncname,'a2x3h_Faxa_rainc','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Faxa_rainc','units','kg m-2 s-1')
        ncwriteatt(ncname,'a2x3h_Faxa_rainc','long_name','Convective precipitation rate')
        ncwriteatt(ncname,'a2x3h_Faxa_rainc','standard_name','convective_precipitation_flux')
        ncwriteatt(ncname,'a2x3h_Faxa_rainc','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Faxa_rainc','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Faxa_rainl','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Faxa_rainl',rainl_assign)
        ncwriteatt(ncname,'a2x3h_Faxa_rainl','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Faxa_rainl','units','kg m-2 s-1')
        ncwriteatt(ncname,'a2x3h_Faxa_rainl','long_name','Large-scale (stable) precipitation rate')
        ncwriteatt(ncname,'a2x3h_Faxa_rainl','standard_name','large_scale_precipitation_flux')
        ncwriteatt(ncname,'a2x3h_Faxa_rainl','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Faxa_rainl','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Faxa_snowc','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Faxa_snowc',snowc_assign)
        ncwriteatt(ncname,'a2x3h_Faxa_snowc','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Faxa_snowc','units','kg m-2 s-1')
        ncwriteatt(ncname,'a2x3h_Faxa_snowc','long_name','Convective snow rate (water equivalent)')
        ncwriteatt(ncname,'a2x3h_Faxa_snowc','standard_name','convective_snowfall_flux')
        ncwriteatt(ncname,'a2x3h_Faxa_snowc','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Faxa_snowc','cell_methods','time: mean')

        nccreate(ncname,'a2x3h_Faxa_snowl','Dimensions',{'a2x3h_nx',1,'a2x3h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x3h_Faxa_snowl',snowl_assign)
        ncwriteatt(ncname,'a2x3h_Faxa_snowl','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x3h_Faxa_snowl','units','kg m-2 s-1')
        ncwriteatt(ncname,'a2x3h_Faxa_snowl','long_name','Large-scale (stable) snow rate (water equivalent)')
        ncwriteatt(ncname,'a2x3h_Faxa_snowl','standard_name','large_scale_snowfall_flux')
        ncwriteatt(ncname,'a2x3h_Faxa_snowl','internal_dname','a2x3h')
        ncwriteatt(ncname,'a2x3h_Faxa_snowl','cell_methods','time: mean')

        ncwriteatt(ncname,'/','creation_method', 'Scaling CAM6 regional product to tower site level met')
        ncwriteatt(ncname,'/','creation_date', datestr(now))
        ncwriteatt(ncname,'/','author', 'Brett Raczka, bmraczka@ucar.edu') 
 

        % Write h1 (1 hour) file
        % 1hr state: a2x1h_Sa_u,a2x1h_Sa_v 
        ncname=[path_scaled_CAM_forcing  enstr{jj} '/CAM6_NR1.cpl_' enstr{jj} '.ha2x1h.' yearstr{ii} '.nc'];
        ncid=netcdf.create(ncname,'CLOBBER');
        netcdf.close(ncid);
        nccreate(ncname,'time','Dimensions',{'time',Inf},'Format','classic');
        ncwrite(ncname,'time',time_1hr)

        if str2num(yearstr{ii})==2011
        ncwriteatt(ncname,'time','units',['days since 2010-07-01 00:00:00'])
        else
        ncwriteatt(ncname,'time','units',['days since 2011-01-01 00:00:00'])
        end
        
        ncwriteatt(ncname,'time','calendar','gregorian')
        ncwriteatt(ncname,'time','bounds','time_bnds')
        ncwriteatt(ncname,'time','cell_methods','time: mean')

        nccreate(ncname,'time_bnds','Dimensions',{'ntb',2,'time',Inf},'Format','classic');
        ncwrite(ncname,'time_bnds',time_bnds_1hr)
        ncwriteatt(ncname,'time_bnds','cell_methods','time: mean')

        nccreate(ncname,'doma_lat','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lat',SITE_lat)
        ncwriteatt(ncname,'doma_lat','_FillValue', 0.)
        ncwriteatt(ncname,'doma_lat','units','degrees north')
        ncwriteatt(ncname,'doma_lat','long_name','latitude')
        ncwriteatt(ncname,'doma_lat','standard_name','latitude')
        ncwriteatt(ncname,'doma_lat','internal_dname','doma')

        nccreate(ncname,'doma_lon','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lon',SITE_lon)
        ncwriteatt(ncname,'doma_lon','_FillValue', 0.)
        ncwriteatt(ncname,'doma_lon','units','degrees east')
        ncwriteatt(ncname,'doma_lon','long_name','longitude')
        ncwriteatt(ncname,'doma_lon','standard_name','longitude')
        ncwriteatt(ncname,'doma_lon','internal_dname','doma')

        nccreate(ncname,'doma_area','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_area',SITE_doma_area)
        ncwriteatt(ncname,'doma_area','_FillValue', 0.)
        ncwriteatt(ncname,'doma_area','units','radian^2')
        ncwriteatt(ncname,'doma_area','long_name','cell_area_model')
        ncwriteatt(ncname,'doma_area','standard_name','cell area from model')
        ncwriteatt(ncname,'doma_area','internal_dname','doma')

        nccreate(ncname,'a2x1h_Sa_u','Dimensions',{'a2x1h_nx',1,'a2x1h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x1h_Sa_u',u_assign)
        ncwriteatt(ncname,'a2x1h_Sa_u','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x1h_Sa_u','units','m s-1')
        ncwriteatt(ncname,'a2x1h_Sa_u','long_name','Zonal wind at the lowest model level')
        ncwriteatt(ncname,'a2x1h_Sa_u','standard_name','eastward_wind')
        ncwriteatt(ncname,'a2x1h_Sa_u','internal_dname','a2x1h')
        ncwriteatt(ncname,'a2x1h_Sa_u','cell_methods','time: mean')

        nccreate(ncname,'doma_aream','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_aream',SITE_doma_aream)
        ncwriteatt(ncname,'doma_aream','_FillValue', 0.)
        ncwriteatt(ncname,'doma_aream','units','radian^2')
        ncwriteatt(ncname,'doma_aream','long_name','cell_area_mapping')
        ncwriteatt(ncname,'doma_aream','standard_name','cell area from mapping file')
        ncwriteatt(ncname,'doma_aream','internal_dname','doma')

        nccreate(ncname,'doma_mask','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_mask',SITE_doma_mask)
        ncwriteatt(ncname,'doma_mask','_FillValue', 0.)
        ncwriteatt(ncname,'doma_mask','units','1')
        ncwriteatt(ncname,'doma_mask','long_name','mask')
        ncwriteatt(ncname,'doma_mask','standard_name','mask')
        ncwriteatt(ncname,'doma_mask','internal_dname','doma')

        nccreate(ncname,'doma_frac','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_frac',SITE_doma_frac)
        ncwriteatt(ncname,'doma_frac','_FillValue', 0.)
        ncwriteatt(ncname,'doma_frac','units','1')
        ncwriteatt(ncname,'doma_frac','long_name','area_fraction')
        ncwriteatt(ncname,'doma_frac','standard_name','area fraction')
        ncwriteatt(ncname,'doma_frac','internal_dname','doma')


        nccreate(ncname,'a2x1h_Sa_v','Dimensions',{'a2x1h_nx',1,'a2x1h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x1h_Sa_v',v_assign)
        ncwriteatt(ncname,'a2x1h_Sa_v','_FillValue', single(1.e30))
        ncwriteatt(ncname,'a2x1h_Sa_v','units','m s-1')
        ncwriteatt(ncname,'a2x1h_Sa_v','long_name','Meridional wind at the lowest model level')
        ncwriteatt(ncname,'a2x1h_Sa_v','standard_name','northward_wind')
        ncwriteatt(ncname,'a2x1h_Sa_v','internal_dname','a2x1h')
        ncwriteatt(ncname,'a2x1h_Sa_v','cell_methods','time: mean')
        ncwriteatt(ncname,'/','creation_method', 'Scaling CAM6 regional product to tower site level met')
        ncwriteatt(ncname,'/','creation_date', datestr(now))
        ncwriteatt(ncname,'/','author', 'Brett Raczka, bmraczka@ucar.edu')
       
      
       end %% ensemble assignment loop  

% Important to clear the CAM6 file variables because leap year makes file size inconsistent

clear Faxa_swndr Faxa_swvdr Faxa_swndf Faxa_swvdf Faxa_rainc Faxa_rainl Faxa_snowc Faxa_snowl 
clear  Sa_u Sa_v Sa_tbot Sa_shum Sa_pbot Faxa_lwdn Sa_z Sa_ptem Sa_dens Sa_pslv Sa_topo 

end   %% Year Loop

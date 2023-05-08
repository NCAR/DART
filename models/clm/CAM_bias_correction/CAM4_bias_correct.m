clear all
close all

% This script applies a bias correction to CAM4 reanalysis meteorology roughly
% following Yong-Fei Zhang (UT-Austin) dissertation approach. Meteorology data
% from the site location of interest should be used as 'truth'.
% All met variables use 'scaling' approach with exception for snow/rain precip.
% *Warning* This script is intended to be used after running CAM4_site_grid.sh

%%%% This section requires user input %%%%%%
% Site location. US-NR1 flux tower used as example
SITE_lat=40.03;
SITE_lon=-105.55+360; % lon should be positive (degrees East only)

%% SITE_doma_area=0.00109327562271889; % m^2 cell area (unique to grid location)
%% SITE_doma_mask=1   ; % all values =1

% Input data files
path_CAM = '/glade/work/bmraczka/CAM4_NR1/';
path_towermet = '/glade/work/bmraczka/SIFMIP2/tower_met_forcing/';
ens_mem=80;  % CAM4 reanalysis provides 80 total members, 1-80 is valid

% Output data files
path_scaled_CAM = '<enter output file path here>';
% Include Diagnostics?  true/false
Diagnostics=true;

% End of user input section
%-----------------------------------------------------------------------------

% Site level extract CAM4:
% One file includes all met forcing.
% noleap calendar
% SOLAR: a2x6h_Faxa_swvdf, a2x6h_Faxa_swndf, a2x6h_Faxa_swvdr, a2x6h_Faxa_swndr, a2x6h_Faxa_lwdn
% PRECIP: a2x6h_Faxa_rainc, a2x6h_Faxa_rainl, a2x6h_Faxa_snowc, a2x6h_Faxa_snowl
% MET:    a2x6h_Sa_pbot, a2x6h_Sa_shum, a2x6h_Sa_tbot, a2x6h_Sa_v, a2x6h_Sa_u
% unchanged: a2x6h_Sa_pslv, a2x6h_Sa_dens, a2x6h_Sa_ptem, a2x6h_Sa_z


% 'yearstr' is list of years where both CAM4 reanalysis and tower met forcing
% is available and the years we will create the bias-correct/scaled reanalysis product
yearstr={'2001','2002','2003','2004','2005','2006','2007','2008','2009','2010'};
% 'yeartower' adds previous year to cover the UTC to MST shift
yeartower={'2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010'};

% Generate ensemble cell array
ens_range=[1:ens_mem];
enstr=cell(1,ens_mem);

for i=1:ens_mem
    enstr{i}=sprintf('%04d', ens_range(i));
end

% Check for existence of output file from CAM4_site_grid.sh
% and set some fixed variables from grid cell collocated with site
site_grid_file = [path_CAM enstr{1} '/CAM4_NR1.cpl_' enstr{1} '.ha2x1dx6h.' yearstr{1} '.nc'];

if exist(site_grid_file,'file')
    SITE_doma_area=ncread(site_grid_file,'doma_area');
    SITE_doma_mask=ncread(site_grid_file,'doma_mask');
else
    error('ERROR !! missing site_grid_file, make sure CAM4_site_grid.sh was run:  %s.',site_grid_file)
end

% Site level met forcing using PLUMBER2 protocol
% time -->30 min increments in LST (MST)
% Tair-->Kelvin
% Qair--> specific humidity (kg kg-1)
% Wind--> (m/s)
% SWdown --> (W/m2)
% LWdown ---> (W/m2)
% Precip --> (mm/s) or (kg m-2 s-1)
% Psurf  --> (Pa)

% Load Tower Met Forcing Data

TBOT_main=load_towermet('Tair',yeartower,path_towermet);
SH_main=load_towermet('Qair',yeartower,path_towermet);
WIND_main=load_towermet('Wind',yeartower,path_towermet);
FSDS_main=load_towermet('SWdown',yeartower,path_towermet);
FLDS_main=load_towermet('LWdown',yeartower,path_towermet);
PRECTmms_main=load_towermet('Precip',yeartower,path_towermet);
PSRF_main=load_towermet('Psurf',yeartower,path_towermet);
YEAR_main=load_towermet('year',yeartower,path_towermet);


% Main Loop: Loads CAM, applies correction, writes corrected CAM file
for ii = 1:length(yearstr);
    
    % Time Zones:
    % Tower Meteorology Forcing (PLUMBER):  LST (MST):  UTC-7
    % CAM4  Forcing: UTC
    
    % CAM4 is 'noleap' calendar whereas tower forcing includes leap days
    % Code accounts for this
    
    % Syncing time steps:
    % Tower Forcing (LST) must be advanced 7 hours to synchronize with CAM4 UTC
    % Total shift of Tower Forcing: +7 hours
    
    
    % Select only indices of YEAR ii
    indices=find(YEAR_main==str2num(yearstr{ii}));
    
    % Pushing met forcing 7 hours forward (go back 14 indices)
    ta_30=TBOT_main(indices(1)-14:indices(end)-14)';                         % TBOT (Kelvin)
    q_30=SH_main(indices(1)-14:indices(end)-14)';                            % SH (kg kg-1)
    wind_30=WIND_main(indices(1)-14:indices(end)-14)';                       % Total Wind (m/s)
    sw_30=FSDS_main(indices(1)-14:indices(end)-14)';                         % FSDS (W/m2)
    lw_30=FLDS_main(indices(1)-14:indices(end)-14)';                         % FLDS (W/m2)
    ppt_30=PRECTmms_main(indices(1)-14:indices(end)-14)';                    % PRECTmms  (mm/s) or (kg m-2 s-1)
    ps_30=PSRF_main(indices(1)-14:indices(end)-14)';                         % PSRF (Pa)
    
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
    
    %% Generate 6 hourly averages to accomodate CAM4
    
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
    
    ppt_3hr=NaN; ta_3hr=NaN; %Placeholders
    clear ta_30 q_30 wind_30 sw_30 lw_30 ppt_30 ps_30
    
    % Load CAM Ensemble Loop
    for jj = 1:ens_mem;
        
        %Downloading all ensemble members for all  reanalysis variables
        
        Faxa_swndr(jj,:) = load_CAM('a2x6h_Faxa_swndr',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Faxa_swvdr(jj,:) = load_CAM('a2x6h_Faxa_swvdr',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Faxa_swndf(jj,:) = load_CAM('a2x6h_Faxa_swndf',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Faxa_swvdf(jj,:) = load_CAM('a2x6h_Faxa_swvdf',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Faxa_rainc(jj,:) = load_CAM('a2x6h_Faxa_rainc',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Faxa_rainl(jj,:) = load_CAM('a2x6h_Faxa_rainl',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Faxa_snowc(jj,:) = load_CAM('a2x6h_Faxa_snowc',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Faxa_snowl(jj,:) = load_CAM('a2x6h_Faxa_snowl',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Sa_u(jj,:) = load_CAM('a2x6h_Sa_u',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Sa_v(jj,:) = load_CAM('a2x6h_Sa_v',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Sa_tbot(jj,:) = load_CAM('a2x6h_Sa_tbot',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Sa_shum(jj,:) = load_CAM('a2x6h_Sa_shum',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Sa_pbot(jj,:) = load_CAM('a2x6h_Sa_pbot',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Faxa_lwdn(jj,:) = load_CAM('a2x6h_Faxa_lwdn',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Sa_z(jj,:) = load_CAM('a2x6h_Sa_z',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Sa_ptem(jj,:) = load_CAM('a2x6h_Sa_ptem',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Sa_dens(jj,:) = load_CAM('a2x6h_Sa_dens',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        Sa_pslv(jj,:) = load_CAM('a2x6h_Sa_pslv',enstr{jj},yearstr{ii},path_CAM,'CAM4');
        
    end %% Load CAM ensemble loop
    
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
    
    % Need this initialization for precip loop
    Faxa_rain=Faxa_rainl;
    Faxa_snow=Faxa_snowl;
    
    
    %% Calculate the scaled correction factors for each year
    %% SOLAR 6 hour resolution
    
    for ind=1:(length(sw_6hr));
        rad_sum=Faxa_swndr_mean(ind)+Faxa_swvdr_mean(ind)+Faxa_swndf_mean(ind)+Faxa_swvdf_mean(ind);
        rad_sum_ratio=(sw_6hr(ind)./(rad_sum));
        if (rad_sum>0)
            Faxa_swndr_scale(ind)= Faxa_swndr_mean(ind) .* rad_sum_ratio;
            Faxa_swvdr_scale(ind)= Faxa_swvdr_mean(ind) .* rad_sum_ratio;
            Faxa_swndf_scale(ind)= Faxa_swndf_mean(ind) .* rad_sum_ratio;
            Faxa_swvdf_scale(ind)= Faxa_swvdf_mean(ind) .* rad_sum_ratio;
        end
        clear rad_sum rad_sum_ratio
    end
    
    %% PRECIP 6  hour resolution
    for ind=1:(length(ppt_6hr));
        
        % Precip loop needs to query adjusted temperature to re-assign to snow or rain
        
        if ta_6hr(ind)>273.15
            % Forcing all unadjusted precip ensemble members to new variable Faxa_rain according to adjusted temperature
            Faxa_rain(:,ind)=Faxa_rainc(:,ind)+Faxa_rainl(:,ind)+Faxa_snowc(:,ind)+Faxa_snowl(:,ind);
            Faxa_snow(:,ind)=0;
            
            % Adjustment will only be applied to total rain. CLM does not care if precip is convective/large scale
            Faxa_rainc_scale(ind)= 0;
            Faxa_rainl_scale(ind)= 0;
            Faxa_snowc_scale(ind)= 0;
            Faxa_snowl_scale(ind)= 0;
            Faxa_rain_scale(ind)= ppt_6hr(ind);
        else
            % Forcing all unadjusted precip ensemble members to new variable Faxa_ snow according to adjusted temperature
            Faxa_snow(:,ind)= Faxa_rainc(:,ind)+Faxa_rainl(:,ind)+Faxa_snowc(:,ind)+Faxa_snowl(:,ind);
            Faxa_rain(:,ind)=0;
            
            % Adjustment will only be applied to total snow. CLM does not care if precip is convective/large scale
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
    %% Put clamping (hard bounds) on non-physical values.
    
    % SOLAR
    % e.g Faxa_swndr(ensemble(ens_mem),time(6hr))
    Faxa_swndr_adjust  = Faxa_swndr+repmat(Faxa_swndr_scale-Faxa_swndr_mean,ens_mem,1); Faxa_swndr_adjust(Faxa_swndr_adjust<0)=0;
    Faxa_swvdr_adjust  = Faxa_swvdr+repmat(Faxa_swvdr_scale-Faxa_swvdr_mean,ens_mem,1); Faxa_swvdr_adjust(Faxa_swvdr_adjust<0)=0;
    Faxa_swndf_adjust  = Faxa_swndf+repmat(Faxa_swndf_scale-Faxa_swndf_mean,ens_mem,1); Faxa_swndf_adjust(Faxa_swndf_adjust<0)=0;
    Faxa_swvdf_adjust  = Faxa_swvdf+repmat(Faxa_swvdf_scale-Faxa_swvdf_mean,ens_mem,1); Faxa_swvdf_adjust(Faxa_swvdf_adjust<0)=0;
    
    %PRECIP
    % Purposely setting convective rain and snow to zero
    % Purposely setting snow and rain adjusted variables to large scale snow (snowl) and rain (rainl)
    % These will carry the ensemble spread for snow/rain
    
    Faxa_rainc_adjust = Faxa_rainc+repmat(Faxa_rainc_scale-Faxa_rainc_mean,ens_mem,1); Faxa_rainc_adjust(:,:)=0;
    Faxa_rainl_adjust = Faxa_rain+repmat(Faxa_rain_scale-Faxa_rain_mean,ens_mem,1); Faxa_rainl_adjust(Faxa_rainl_adjust<0)=0;
    Faxa_snowc_adjust = Faxa_snowc+repmat(Faxa_snowc_scale-Faxa_snowc_mean,ens_mem,1); Faxa_snowc_adjust(:,:)=0;
    Faxa_snowl_adjust = Faxa_snow+repmat(Faxa_snow_scale-Faxa_snow_mean,ens_mem,1); Faxa_snowl_adjust(Faxa_snowl_adjust<0)=0;
    
    % meridional and zonal winds allowed to go negative
    Sa_u_adjust = Sa_u+repmat(Sa_u_scale-Sa_u_mean,ens_mem,1);
    Sa_v_adjust = Sa_v+repmat(Sa_v_scale-Sa_v_mean,ens_mem,1);
    
    % Remaining variables
    Sa_tbot_adjust = Sa_tbot+repmat(Sa_tbot_scale-Sa_tbot_mean,ens_mem,1); Sa_tbot_adjust(Sa_tbot_adjust<0)=0;
    Sa_shum_adjust = Sa_shum+repmat(Sa_shum_scale-Sa_shum_mean,ens_mem,1); Sa_shum_adjust(Sa_shum_adjust<0)=0;
    Sa_pbot_adjust = Sa_pbot+repmat(Sa_pbot_scale-Sa_pbot_mean,ens_mem,1); Sa_pbot_adjust(Sa_pbot_adjust<0)=0;
    Faxa_lwdn_adjust = Faxa_lwdn+repmat(Faxa_lwdn_scale-Faxa_lwdn_mean,ens_mem,1); Faxa_lwdn_adjust(Faxa_lwdn_adjust<0)=0;
    
    % (3hr STATE Unchanged)
    %Sa_z = Sa_z+repmat(Sa_z_scale-Sa_z_mean,ens_mem,1);
    %Sa_ptem = Sa_ptem+repmat(Sa_ptem_scale-Sa_ptem_mean,ens_mem,1);
    %Sa_dens = Sa_dens+repmat(Sa_dens_scale-Sa_dens_mean,ens_mem,1);
    %Sa_pslv = Sa_pslv+repmat(Sa_pslv_scale-Sa_pslv_mean,ens_mem,1);
    %Sa_topo = Sa_topo+repmat(Sa_topo_scale-Sa_topo_mean,ens_mem,1);
    
    
    if Diagnostics==true
        
        plot_diagnostic(yearstr{ii},Faxa_swndr,Faxa_swndf,Faxa_swvdr,Faxa_swvdf, ...
            Faxa_rainl,Faxa_rainc,Faxa_snowl,Faxa_snowc,Sa_tbot, ...
            Sa_shum,Sa_u,Sa_v,Sa_pbot,Faxa_lwdn,Faxa_rain,Faxa_snow, ...
            Faxa_swndr_adjust,Faxa_swndf_adjust,Faxa_swvdr_adjust, ...
            Faxa_swvdf_adjust,Faxa_rainl_adjust,Faxa_snowl_adjust, ...
            Sa_tbot_adjust,Sa_shum_adjust,Sa_u_adjust,Sa_v_adjust, ...
            Sa_pbot_adjust,Faxa_lwdn_adjust,'CAM4', ...
            sw_1hr,ppt_1hr,ta_1hr,q_1hr,wind_1hr,ps_1hr,lw_1hr, ...
            ppt_3hr,ta_3hr,wind_6hr,ta_6hr,ppt_6hr,sw_6hr);
    end
    
    
    % Assign the adjusted CAM variables to netcdf files for each ensemble member/ year
    
    time_var= ncread([path_CAM enstr{1} '/CAM4_NR1.cpl_' enstr{1} '.ha2x1dx6h.' yearstr{ii} '.nc'],'time');
    dim_time=length(time_var);        % ha2x1dx6h
    time_bnds = ncread([path_CAM enstr{1} '/CAM4_NR1.cpl_' enstr{1} '.ha2x1dx6h.' yearstr{ii} '.nc'],'time_bnds');
    
    
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
    
    for jj = 1:ens_mem;  %% Ensemble assignment loop
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
        
        ncname=[path_scaled_CAM 'CAM4_NR1.cpl_' enstr{jj} '.ha2x1dx6h.' yearstr{ii} '.nc'];
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
        write_netcdf_att(ncname,'doma_area',0.,'m^2','undefined','cell area','doma')
        
        nccreate(ncname,'doma_lat','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lat',SITE_lat)
        write_netcdf_att(ncname,'doma_lat',0.,'degrees north','undefined','latitude','doma')
        
        nccreate(ncname,'doma_lon','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_lon',SITE_lon)
        write_netcdf_att(ncname,'doma_lon',0.,'degrees east','undefined','longitude','doma')
        
        nccreate(ncname,'doma_mask','Dimensions',{'doma_nx',1,'doma_ny',1},'Format','classic');
        ncwrite(ncname,'doma_mask',SITE_doma_mask)
        write_netcdf_att(ncname,'doma_mask',0.,'unitless','undefined','mask','doma')
        
        nccreate(ncname,'a2x6h_Faxa_swndf','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_swndf',swndf_assign)
        write_netcdf_att(ncname,'a2x6h_Faxa_swndf',single(1.e30),'W m-2','Diffuse near-infrared incident solar radiation', ...
            'surface_downward_diffuse_shortwave_flux_due_to_near_infrared_radiation','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Faxa_swndr','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_swndr',swndr_assign)
        write_netcdf_att(ncname,'a2x6h_Faxa_swndr',single(1.e30),'W m-2','Direct near-infrared incident solar radiation', ...
            'surface_downward_direct_shortwave_flux_due_to_near_infrared_radiation','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Faxa_swvdf','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_swvdf',swvdf_assign)
        write_netcdf_att(ncname,'a2x6h_Faxa_swvdf',single(1.e30),'W m-2','Diffuse visible incident solar radiation', ...
            'surface_downward_diffuse_shortwave_flux_due_to_visible_radiation','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Faxa_swvdr','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_swvdr',swvdr_assign)
        write_netcdf_att(ncname,'a2x6h_Faxa_swvdr',single(1.e30),'W m-2','Direct visible incident solar radiation', ...
            'surface_downward_direct_shortwave_flux_due_to_visible_radiation','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_z','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_z',z_assign)
        write_netcdf_att(ncname,'a2x6h_Sa_z',single(1.e30),'m','Height at the lowest model level','height','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_tbot','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_tbot',tbot_assign)
        write_netcdf_att(ncname,'a2x6h_Sa_tbot',single(1.e30),'K','Temperature at the lowest model level','air_temperature', ...
            'a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_ptem','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_ptem',ptem_assign)
        write_netcdf_att(ncname,'a2x6h_Sa_ptem',single(1.e30),'K','Potential temperature at the lowest model level', ...
            'air_potential_temperature','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_shum','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_shum',shum_assign)
        write_netcdf_att(ncname,'a2x6h_Sa_shum',single(1.e30),'kg kg-1','Specific humidity at the lowest model level', ...
            'specific humidity','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_dens','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_dens',dens_assign)
        write_netcdf_att(ncname,'a2x6h_Sa_dens',single(1.e30),'kg m-3','Air density at the lowest model level', ...
            'air_density','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_pbot','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_pbot',pbot_assign)
        write_netcdf_att(ncname,'a2x6h_Sa_pbot',single(1.e30),'Pa','Pressure at the lowest model level', ...
            'air_pressure','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_pslv','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_pslv',pslv_assign)
        write_netcdf_att(ncname,'a2x6h_Sa_pslv',single(1.e30),'Pa','Sea level pressure','air_pressure_at_sea_level', ...
            'a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Faxa_lwdn','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_lwdn',lwdn_assign)
        write_netcdf_att(ncname,'a2x6h_Faxa_lwdn',single(1.e30),'W m-2','Downward longwave heat flux', ...
            'downwelling_longwave_flux','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Faxa_rainc','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_rainc',rainc_assign)
        write_netcdf_att(ncname,'a2x6h_Faxa_rainc',single(1.e30),'kg m-2 s-1','Convective precipitation rate', ...
            'convective_precipitation_flux','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Faxa_rainl','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_rainl',rainl_assign)
        write_netcdf_att(ncname,'a2x6h_Faxa_rainl',single(1.e30),'kg m-2 s-1','Large-scale (stable) precipitation rate', ...
            'large_scale_precipitation_flux','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Faxa_snowc','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_snowc',snowc_assign)
        write_netcdf_att(ncname,'a2x6h_Faxa_snowc',single(1.e30),'kg m-2 s-1','Convective snow rate (water equivalent)', ...
            'convective_snowfall_flux','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Faxa_snowl','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Faxa_snowl',snowl_assign)
        write_netcdf_att(ncname,'a2x6h_Faxa_snowl',single(1.e30),'kg m-2 s-1','Large-scale (stable) snow rate (water equivalent)', ...
            'large_scale_snowfall_flux','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_u','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_u',u_assign)
        write_netcdf_att(ncname,'a2x6h_Sa_u',single(1.e30),'m s-1','Zonal wind at the lowest model level', ...
            'eastward_wind','a2x6h','time: mean')
        
        nccreate(ncname,'a2x6h_Sa_v','Dimensions',{'a2x6h_nx',1,'a2x6h_ny',1,'time',Inf},'Datatype','single','Format','classic');
        ncwrite(ncname,'a2x6h_Sa_v',v_assign)
        write_netcdf_att(ncname,'a2x6h_Sa_v',single(1.e30),'m s-1','Meridional wind at the lowest model level', ...
            'northward_wind','a2x6h','time: mean')
        
        ncwriteatt(ncname,'/','creation_method', 'Scaling CAM4 regional product (ds199.1) to tower site level met')
        ncwriteatt(ncname,'/','creation_date', datestr(now))
        ncwriteatt(ncname,'/','author', 'Brett Raczka, bmraczka@ucar.edu')
        
    end %% Ensemble assignment loop
    
    % Important to clear the CAM4 file variables because leap year makes file size inconsistent
    
    clear Faxa_swndr Faxa_swvdr Faxa_swndf Faxa_swvdf Faxa_rainc Faxa_rainl Faxa_snowc Faxa_snowl
    clear  Sa_u Sa_v Sa_tbot Sa_shum Sa_pbot Faxa_lwdn Sa_z Sa_ptem Sa_dens Sa_pslv
    
end   %% Main Loop

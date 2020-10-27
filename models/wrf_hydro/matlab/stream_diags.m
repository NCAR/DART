
%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id: stream_diags.m $

clear 
close all
clc

diag_dir = 'bucket2/';

%% DATA
nc_state_mean    = strcat(diag_dir,'all_preassim_mean.nc');
nc_state_sd      = strcat(diag_dir,'all_preassim_sd.nc');
nc_inflate_mean  = strcat(diag_dir,'all_output_priorinf_mean.nc');
nc_inflate_sd    = strcat(diag_dir,'all_output_priorinf_sd.nc');
nc_routelink     = strcat(diag_dir,'RouteLink.nc');

nc_obs_output    = strcat(diag_dir,'obs_diag_output.nc');
obs_epoc_file    = strcat(diag_dir,'obs_epoch_001.nc');

% Retrieve the dimensions from the netcdf file
ncid = netcdf.open(nc_obs_output,'NC_NOWRITE');
[~, rbins]  = netcdf.inqDim(ncid, 14); ens_size = rbins-1;
netcdf.close(ncid);

ncid = netcdf.open(nc_state_mean,'NC_NOWRITE');
[~, da_cycles] = netcdf.inqDim(ncid, 0);
[~, n_links]   = netcdf.inqDim(ncid, 1);
LINKS          = 1:n_links;
netcdf.close(ncid);


%% TIME HANDLING:
Time        = ncread(nc_state_mean, 'time'); Nt = length(Time);
origin      = datenum(1601, 1, 1);
current     = Time + origin; 
goodtime    = datestr(current, 'mmm dd'); 
dateformat  = 'dd-mmm-yyyy HH:MM';
detailtime  = datestr(current, dateformat); 

O_time      = double(ncread(obs_epoc_file, 'time')) + origin;

time_label  = ceil( [1, Nt/4, Nt/2, 3*Nt/4, Nt] );
xticks      = current(time_label);
xtickslabel = goodtime(time_label, :);

% TODO: Need a better strategy to handle hurrican events: 
% Interesting events
if n_links < 1000 && da_cycles > 2000
    % SIXMILE CASE:
    % *************
    events = [250, 1660]; % hard-coded: can't find info online!
    
elseif n_links > 1000 && da_cycles > 800
    % FLORENCE CASE:
    % **************
        % Total fatalities: 53
        % Highest wind speed: 137 mph
        % Dates: Aug 31, 2018 ? Sep 19, 2018
        % Date: August 31, 2018 ? September 19, 2018
        % Category: Category 4 major hurricane (NHC/CPHC)
        % Affected areas: North Carolina, South Carolina, Maryland
    
    landfall = '11-Sep-2018 08:00';
    landfall_date = datenum(landfall, dateformat);
    
    events = [1, find(current == landfall_date)];
    
else
    % not so much interesting, can't do nothing ..
    events = [1, size(avg_all_links, 2)];
end


%% PROCESSING DATA:
gauge_id = strtrim(ncread(nc_routelink, 'gages')');

special_gauges_ids  = [02071530, 02086849, 02070500, 0209553650, 02088383, 02047000, 02082585, 02083500 ] ;
special_gauges_cell = {'02071530', '02086849', '02070500', '0209553650', ...
                       '02088383', '02047000', '02082585', '02083500'};

Lg = length(special_gauges_ids); 
k = 0;
for i = 1:n_links
    ob_id = gauge_id(i, :);
    
    if sum( isspace(ob_id) ) ~= 10 
        for t = 1:Lg 
            if str2double(ob_id) == special_gauges_ids(t)
                k = k +1;
                keepInd(k) = i;
                
            end
        end
    end
    
end

gauge_name    = strtrim(ncread(nc_obs_output, 'ObservationTypes')');
char_indices  = gauge_name(:, 13:end);
gauge_num     = length(char_indices);
gauge_index   = zeros(1, gauge_num);
for l = 1:gauge_num
    gauge_index(l) = str2double(char_indices(l,:));
end

for t = 1:Lg
    toPlot(t) = find(gauge_index == keepInd(t));
end

% Reading:
% i. obs_diag output 
  
% ii. obs_seq_to_netcdf output
ensemble = double(ncread(obs_epoc_file, 'observations'));
obs_ind  = -1 * double(ncread(obs_epoc_file, 'obs_type'));

forecast    = cell(8, gauge_num);
observation = cell(7, gauge_num);
analysis    = cell(7, gauge_num);

for i = 1:gauge_num

    k = gauge_index(i);
    
    find_obs    = k == obs_ind;
    
    tmp.obs_val    = ensemble(1,         find_obs);
    tmp.ens_mean_f = ensemble(2,         find_obs);
    tmp.ens_sd_f   = ensemble(4,         find_obs);
    tmp.ensemble_f = ensemble(6:2:end-1, find_obs);
    tmp.obs_var    = ensemble(end,       find_obs); 
    tmp.ens_mean_a = ensemble(3,         find_obs);
    tmp.ens_sd_a   = ensemble(5,         find_obs);
    tmp.ensemble_a = ensemble(7:2:end,   find_obs);
    
    ens_time   = zeros(1, Nt+1);
    Found_time = O_time(find_obs); 
    for j = 1:Nt
        get_t_index = sum(find(current(j) == Found_time));
        
        if get_t_index > 0 % Found data for obs-type k at time j
            ens_time(j) = get_t_index - 1;
        else % No gauge data for this obs-type k at time j
            ens_time(j) = NaN;
        end
    end
    ens_time(   1) = 1; 
    ens_time(Nt+1) = length(Found_time);
    
    obs_val    = zeros(1, Nt);
    ens_mean_f = zeros(1, Nt);
    ens_sd_f   = zeros(1, Nt);
    ensemble_f = zeros(ens_size, Nt);
    obs_var    = zeros(1, Nt);
    ens_mean_a = zeros(1, Nt);
    ens_sd_a   = zeros(1, Nt);
    ensemble_a = zeros(ens_size, Nt);
    for j = 1:Nt
        if ~isnan(ens_time(j)) && ~isnan(ens_time(j+1))
            obs_val(j)       = mean(tmp.obs_val   (:, ens_time(j):ens_time(j+1)), 2);
            ens_mean_f(j)    = mean(tmp.ens_mean_f(:, ens_time(j):ens_time(j+1)), 2);
            ens_sd_f(j)      = mean(tmp.ens_sd_f  (:, ens_time(j):ens_time(j+1)), 2);
            ensemble_f(:, j) = mean(tmp.ensemble_f(:, ens_time(j):ens_time(j+1)), 2);
            obs_var(j)       = mean(tmp.obs_var   (:, ens_time(j):ens_time(j+1)), 2);
            ens_mean_a(j)    = mean(tmp.ens_mean_a(:, ens_time(j):ens_time(j+1)), 2);
            ens_sd_a(j)      = mean(tmp.ens_sd_a  (:, ens_time(j):ens_time(j+1)), 2);
            ensemble_a(:, j) = mean(tmp.ensemble_a(:, ens_time(j):ens_time(j+1)), 2);
            
        elseif ~isnan(ens_time(j))
            obs_val(j)       = tmp.obs_val   (:, ens_time(j));
            ens_mean_f(j)    = tmp.ens_mean_f(:, ens_time(j));
            ens_sd_f(j)      = tmp.ens_sd_f  (:, ens_time(j));
            ensemble_f(:, j) = tmp.ensemble_f(:, ens_time(j));
            obs_var(j)       = tmp.obs_var   (:, ens_time(j));
            ens_mean_a(j)    = tmp.ens_mean_a(:, ens_time(j));
            ens_sd_a(j)      = tmp.ens_sd_a  (:, ens_time(j));
            ensemble_a(:, j) = tmp.ensemble_a(:, ens_time(j));
            
        elseif ~isnan(ens_time(j+1))
            obs_val(j)       = tmp.obs_val   (:, ens_time(j+1));
            ens_mean_f(j)    = tmp.ens_mean_f(:, ens_time(j+1));
            ens_sd_f(j)      = tmp.ens_sd_f  (:, ens_time(j+1));
            ensemble_f(:, j) = tmp.ensemble_f(:, ens_time(j+1));
            obs_var(j)       = tmp.obs_var   (:, ens_time(j+1));
            ens_mean_a(j)    = tmp.ens_mean_a(:, ens_time(j+1));
            ens_sd_a(j)      = tmp.ens_sd_a  (:, ens_time(j+1));
            ensemble_a(:, j) = tmp.ensemble_a(:, ens_time(j+1));
            
        else
            obs_val(j)       = NaN;
            ens_mean_f(j)    = NaN;
            ens_sd_f(j)      = NaN;
            ensemble_f(:, j) = NaN;
            obs_var(j)       = NaN;
            ens_mean_a(j)    = NaN;
            ens_sd_a(j)      = NaN;
            ensemble_a(:, j) = NaN;
            
        end
    end
    clear tmp
        
    varname_f = strcat(gauge_name(i,:), '_guess');
    varname_a = strcat(gauge_name(i,:), '_analy');
    
    % Manage the forecast copies
    tmp_f       = squeeze( double(ncread(nc_obs_output, varname_f)) );
    rmse_f      = tmp_f(7, :);  %abs(ens_mean_f - obs_val);
    bias_f      = tmp_f(8, :);  %ens_mean_f - obs_val;
    totspread_f = tmp_f(10, :); %sqrt( ens_sd_f.^2 + obs_var.^2 );
    rank_hist_f = squeeze( double(ncread(nc_obs_output, strcat(varname_f, '_RankHist'))) );
    
    % Manage the observation copies
    obs_poss    = tmp_f(1, :);
    obs_used    = tmp_f(2, :);
    obs_w_qc7   = tmp_f(22, :);
    
    % Manage the analysis copies
    tmp_a       = squeeze( double(ncread(nc_obs_output, varname_a)) );
    rmse_a      = tmp_a(7, :);  %abs(ens_mean_a - obs_val);
    bias_a      = tmp_a(8, :);  %ens_mean_a - obs_val;
    totspread_a = tmp_a(10, :); %sqrt( ens_sd_a.^2 + obs_var.^2 );
    
    
    % Construct the forecast, observation and analysis cells.
    forecast(:, i)    = { varname_f      ; ...
                          rmse_f         ; ...
                          bias_f         ; ... 
                          ens_sd_f       ; ...
                          totspread_f    ; ...
                          ens_mean_f     ; ...
                          rank_hist_f    ; ...
                          ensemble_f     ; ...
                        };
                 
    observation(:, i) = { k              ; ...
                          gauge_id(k, :) ; ...
                          obs_poss       ; ...
                          obs_used       ; ...
                          obs_w_qc7      ; ...
                          obs_val        ; ...
                          obs_var        ; ...
                        };  
                          
                 
    analysis(:, i)    = { varname_a      ; ...
                          rmse_a         ; ...
                          bias_a         ; ... 
                          ens_sd_a       ; ...
                          totspread_a    ; ...
                          ens_mean_a     ; ...
                          ensemble_a     ; ...
                        };
end


%% READ STATE AND INFLATION IN SPACE:
avg_all_links = ncread(nc_state_mean, 'qlink1');
sd_all_links  = ncread(nc_state_sd  , 'qlink1');

% Get value of the mean and sd of the state
sf_mean = squeeze(avg_all_links(:, events));
sf_sd   = squeeze(sd_all_links(:, events));

% inflation mean and sd:
inf_mean = ncread(nc_inflate_mean, 'qlink1');
inf_sd   = ncread(nc_inflate_sd, 'qlink1'); 


%% MISC:
gR= [ 0,153,0 ]/255; 
bL= [30,144,255]/255;
rD= [ 255,51,51 ]/255;
gY= [ 150,150,150 ]/255; 
pR= [ 153,51,255 ]/255;
oR= [ 255,153,51 ]/255;


%% DISPLAY RESULTS: %
% ***************** %

%% TIME SERIES EVOLUTION:

for t = 1:Lg
    figure('uni','pi','pos',[200, 600, 600, 900])
    
    o = toPlot(t);
    
    subplot(311)
    en = plot(current, forecast{8, o}   , ':', 'Color', gY); hold on
    ob = plot(current, observation{6, o}, '.', 'Color', gR); grid on 
    mf = plot(current, forecast{6, o}   , '-', 'Color', 'k'); 
    ma = plot(current, analysis{6, o}   , '-', 'Color', 'b'); 
    
    set(gca, 'FontSize', 14, 'XTick', xticks, 'XTickLabel', xtickslabel, 'YGrid', 'off')
    ylabel('Stream flow (cms)', 'FontSize', 14)
    
    legend([ob, en(1), mf, ma], 'Observation', 'Prior Members', 'Prior Mean', 'Posterior Mean', ...
                                'Location', 'Best')
                            
    str1 = gauge_name(o,:);
    str2 = 'Time-Series: Obs, Prior Mean and Ensemble';
    title({str1, str2}, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none')
    
    
    subplot(312)
    yyaxis left 

    rf = plot(current, forecast{2, o}, '-', 'Color', 'k'); hold on 
    sf = plot(current, forecast{4, o}, '-', 'Color', 'b'); grid on
    ra = plot(current, analysis{2, o}, '--', 'Color', 'k'); hold on 
    sa = plot(current, analysis{4, o}, '--', 'Color', 'b'); grid on

    set(gca, 'YColor', 'k', 'FontSize', 14, 'XTick', xticks, 'XTickLabel', xtickslabel, 'YGrid', 'off')

    ylabel('Error, Spread (cms)', 'FontSize', 14)

    yyaxis right

    op = plot(current, observation{3, o}, 'o', 'Color', 'r'); hold on
    ou = plot(current, observation{4, o}, 'x', 'Color', 'r');
    set(gca, 'YColor', 'k')
    
    legend([rf, sf, ra, sa, op, ou],    [ 'Prior RMSE, Avg: ' num2str(nanmean(forecast{2, o})) ], ... 
                                        [ 'Prior spread, Avg: '  num2str(nanmean(forecast{4, o})) ], ...
                                        [ 'Posterior RMSE, Avg: ' num2str(nanmean(analysis{2, o})) ], ... 
                                        [ 'Posterior spread, Avg: '  num2str(nanmean(analysis{4, o})) ], ...
                                        'Possible Obs', 'Used Obs', ...
                                        'Location', 'Best')

    ylabel([ 'Number of obs (# rejected: ' num2str(nansum(observation{3, o} - observation{4, o})) ')' ] , 'FontSize', 14)

    str1 = sprintf('Gauge ID: %s', observation{2, o});
    str2 = 'Time-Series: Prior and Posterior Statistics';
    title({str1, str2}, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none')


    subplot(313)
    B = bar(1:ens_size+1, mean(forecast{7, o}, 2)); grid on
    B.FaceColor = gY;
    B.BarWidth  = 1.;

    set(gca, 'FontSize', 12, 'XLim', [1, ens_size+1], 'XTick', (10:10:ens_size-10))
    xlabel('Observation Rank (among ensemble members)', 'FontSize', 14)
    ylabel('Count', 'FontSize', 14)

    str1 = 'Rank Histogram';
    str2 = sprintf('Binning over entire %d DA cycles', da_cycles);
    title({str1, str2}, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none')

end


%% All obs, mean scatter
figure('uni','pi','pos',[200, 400, 600, 500])

prior_mean = zeros(1, gauge_num);
post_mean  = zeros(1, gauge_num);
obs_mean   = zeros(1, gauge_num);
for o = 1:gauge_num
    prior_mean(o) = nanmean(forecast{6, o});
    post_mean(o)  = nanmean(analysis{6, o});
    obs_mean(o)   = nanmean(observation{6, o});
end
sMax = 1.e2;

plot(obs_mean, prior_mean, 'ob'); hold on 
plot(obs_mean, post_mean, 'or'); grid on 

plot([0, sMax], [0, sMax], '-k')

set(gca, 'FontSize', 14, 'XLim', [0, sMax], 'YLim', [0, sMax])
ylabel('Observations (cms)', 'FontSize', 16)
xlabel('Streamflow Ensemble Mean (cms)', 'FontSize', 16)

legend('Prior', 'Posterior')

title('Streamflow Observations vs Filtered Ensemble Means', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none')


%% Taylor
figure('uni','pi','pos',[200, 400, 600, 500])

for t = 1:Lg
    
    obs_series = observation{6, toPlot(t)};
    pri_series = forecast{6, toPlot(t)};
    pos_series = analysis{6, toPlot(t)};
    
    tS{t} = [obs_series', pri_series', pos_series'];
end

data_locs = special_gauges_cell; 
data_tags = {'Prior', 'Posterior'};

norm_taylor_diag(tS, data_locs, data_tags);


%% Number of assimilated and possible obs + Spread
figure('uni','pi','pos',[20, 400, 1800, 400])

all_obs_pos = zeros(gauge_num, 1);
all_obs_use = zeros(gauge_num, 1);
all_pri_aes = zeros(gauge_num, 1);
all_pos_aes = zeros(gauge_num, 1);
for o = 1:gauge_num
    all_obs_pos(o) = nansum(observation{3, o});
    all_obs_use(o) = nansum(observation{4, o});
    
    all_pri_aes(o) = nanmean(forecast{4, o});
    all_pos_aes(o) = nanmean(analysis{4, o});
end
merge_obs_num = [all_obs_use, all_obs_pos] / 1000;

yyaxis left
B_obs = bar(merge_obs_num, 'stacked'); grid on

gaugelab = [1, 22:20:gauge_num]; colormap autumn

set(gca, 'YColor', 'k', 'FontSize', 14, 'XLim', [0, gauge_num+1], ...
         'XTick', gaugelab, 'XTickLabel', char_indices(gaugelab, :))

ylabel('Number of observations x 10^3', 'FontSize', 16)
xlabel('Gauge ID', 'FontSize', 16)

yyaxis right
pri_aes = plot(1:gauge_num, all_pri_aes, 'ob', ...
          'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on 
pos_aes = plot(1:gauge_num, all_pos_aes, 'ok', ...
          'MarkerFaceColor', 'k', 'MarkerSize', 6); 

set(gca, 'YColor', 'k', 'YScale', 'log')
ylabel('Ensemble Spread', 'FontSize', 16)

title('Available Observations and Ensemble Spread at Streamflow Gauges', ...
      'FontSize', 20, 'FontWeight', 'Normal')

legend([B_obs(1), B_obs(2), pri_aes, pos_aes], ...
       [ 'Used obs: ' num2str(sum(all_obs_use)) ], ...
       [ 'Possible obs: ' num2str(sum(all_obs_pos)) ], ...
       'Prior Ensemble Spread', 'Posterior Ensemble Spread', ...
       'Location', 'North')
   

%% STATE IN SPACE:
figure('uni','pi','pos', [10, 300, 800, 1000])

subplot(211)
plot_connections(sf_mean(:, 2), get(gca, 'position'), 0, 'Streamflow Mean (cms)')
title([ 'Event: ' detailtime(events(2), :) ], 'FontSize', 15, 'FontWeight', 'normal')

subplot(212)
plot_connections(sf_sd(:, 2), get(gca, 'position'), 0, 'Streamflow SD (cms)')
title([ 'Event: ' detailtime(events(2), :) ], 'FontSize', 15, 'FontWeight', 'normal')

%% Inflation 
figure('uni','pi','pos', [10, 300, 800, 1000])

subplot(211)
plot_connections(mean(inf_mean, 2), get(gca, 'position'), 0, 'Inflation Mean (cms)')
title('Temperal Average of Inflation Mean', 'FontSize', 15, 'FontWeight', 'normal')

subplot(212)
plot_connections(mean(inf_sd, 2), get(gca, 'position'), 0, 'Inflation SD (cms)')
title('Temperal Average of Inflation Standard Deviation', 'FontSize', 15, 'FontWeight', 'normal')

% <next few lines under version control, do not edit>
% $URL: $
% $Revision: $
% $Date: $

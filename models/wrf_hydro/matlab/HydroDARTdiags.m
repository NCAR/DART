function [observation, openloop, forecast, analysis, exp] = ...
    HydroDARTdiags(dir_exps, obs, dir_ol, disp_res, plot_state, fig_to_pdf)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
% 
% ** dir_exps: Experiment directories 
%
% ** obs: The guages where statistics are displayed.
%         These have to be inside the domain area. 
%         e.g.: [02086849, 0208521324, ...] or simply 0 (will display all
%         available gauges)
%
% ** dir_ol: Optional Open Loop directory
%
% ** disp_res: Toggle to display the results
%
% ** plot_state: Option to display the entire state on the link network
%
% ** output are structures for all time-series data
%
%   example: dir_exps   = {'da_exp1', 'da_exp2', 'da_exp3'};
%            obs        = [02086849, 0208521324, 02085000];
%            dir_ol     = {'openloop_exp'};
%            disp_res   = 1;
%            plot_state = 1;
%
%            HydroDARTdiags(dir_exps, obs, dir_ol, disp_res, plot_state)
%
% DART $Id: HydroDARTdiags.m $

if nargin < 1
    error('No arguments provided!! Please enter the following: (1) Name of experiment directories and (2) Gauges to diagnose')
    
elseif nargin > 2
    plot_ol = true;
    
elseif nargin == 2 % no open-loop case
    plot_ol    = false; 
    plot_state = false;
    disp_res   = 1;
    fig_to_pdf = 'results.pdf';
end

gY = [ 150, 150, 150 ]/255;
lB = [ 153, 255, 255 ]/255;
bK = [   0,   0,   0 ]/255;
bL = [  30, 144, 255 ]/255;
rD = [ 255,  51,  51 ]/255;
gR = [   0, 153,   0 ]/255;
oR = [ 255, 153,  51 ]/255;


%% DATA
num_exps = length(dir_exps);

if num_exps > 4
    warning('Too many experiments produce a lot of "cluttered" figures!!')
end

% prior or posterior inflation 
isprior = '/all_output_priorinf_mean.nc'; 
ispost  = '/all_output_postinf_mean.nc'; 

inf_flavor = cell(num_exps, 2);
inf_flav_n = zeros(num_exps, 2);
for e = 1:num_exps
    if exist(char(strcat(dir_exps(e), isprior)), 'file') == 2 && ...
        exist(char(strcat(dir_exps(e), ispost)), 'file') == 2
        inf_flavor{e, 1} = '/all_output_priorinf_';
        inf_flavor{e, 2} = '/all_output_postinf_';
        inf_flav_n(e, 1) = 1;
        inf_flav_n(e, 2) = 1;
        
    elseif exist(char(strcat(dir_exps(e), isprior)), 'file') == 2
        inf_flavor{e, 1} = '/all_output_priorinf_'; inf_flavor{e, 2} = '';
        inf_flav_n(e, 1) = 1; inf_flav_n(e, 2) = 0;
        
    elseif exist(char(strcat(dir_exps(e), ispost)), 'file') == 2
        inf_flavor{e, 1} = ''; inf_flavor{e, 2} = '/all_output_postinf_';
        inf_flav_n(e, 1) = 0;  inf_flav_n(e, 2) = 1;
    end
end

% Do we have any hybrid weight files?
ishybrid = '/all_output_hybridweight_mean.nc'; 

hyb_flavor = cell(num_exps);
hyb_flav_n = zeros(num_exps);
for e = 1:num_exps
    if exist(char(strcat(dir_exps(e), ishybrid)), 'file') == 2
        hyb_flavor{e} = '/all_output_hybridweight_'; 
        hyb_flav_n(e) = 1; 
    end
end

nc = struct;
for e = 1:num_exps
    
    diag_dir = dir_exps(e);
    
    nc(e).state_mean_pr = char(strcat(diag_dir, '/all_preassim_mean.nc'           )); % aggregated prior state_mean
    nc(e).state_sd_pr   = char(strcat(diag_dir, '/all_preassim_sd.nc'             )); % aggregated prior state_sd
    
    nc(e).state_mean_po = char(strcat(diag_dir, '/all_analysis_mean.nc'           )); % aggregated analysis state_mean
    nc(e).state_sd_po   = char(strcat(diag_dir, '/all_analysis_sd.nc'             )); % aggregated analysis state_sd

    if inf_flav_n(e, 1) > 0 && inf_flav_n(e, 2) > 0
        nc(e).pr_inflate_mean  = char(strcat(diag_dir, inf_flavor{e, 1}, 'mean.nc')); % aggregated inf_mean
        nc(e).pr_inflate_sd    = char(strcat(diag_dir, inf_flavor{e, 1}, 'sd.nc'  )); % aggregated inf_std 
        nc(e).po_inflate_mean  = char(strcat(diag_dir, inf_flavor{e, 2}, 'mean.nc')); % aggregated inf_mean
        nc(e).po_inflate_sd    = char(strcat(diag_dir, inf_flavor{e, 2}, 'sd.nc'  )); % aggregated inf_std 
        
    elseif inf_flav_n(e, 1) > 0
        nc(e).pr_inflate_mean  = char(strcat(diag_dir, inf_flavor{e, 1}, 'mean.nc')); % aggregated inf_mean
        nc(e).pr_inflate_sd    = char(strcat(diag_dir, inf_flavor{e, 1}, 'sd.nc'  )); % aggregated inf_std 
        
    elseif inf_flav_n(e, 2) > 0
        nc(e).po_inflate_mean  = char(strcat(diag_dir, inf_flavor{e, 2}, 'mean.nc')); % aggregated inf_mean
        nc(e).po_inflate_sd    = char(strcat(diag_dir, inf_flavor{e, 2}, 'sd.nc'  )); % aggregated inf_std 
    end
 
    if hyb_flav_n(e) > 0
        nc(e).hybrid_mean = char(strcat(diag_dir, hyb_flavor{e}, 'mean.nc')); % aggregated hyb_mean
        nc(e).hybrid_sd   = char(strcat(diag_dir, hyb_flavor{e}, 'sd.nc'  )); % aggregated hyb_std 
    end

    nc(e).routelink     = char(strcat(diag_dir, '/../RouteLink.nc'                   )); % routelink file

    nc(e).obs_diag      = char(strcat(diag_dir, '/obs_diag_output.nc'             )); % output of obs_diag
    nc(e).obs_epoc      = char(strcat(diag_dir, '/obs_epoch_001.nc'               )); % output of obs_seq_to_netcdf
    
    % open loop
    if plot_ol 
        ol.obs_diag     = char(strcat(dir_ol  , '/obs_diag_output.nc'             )); 
        ol.obs_epoc     = char(strcat(dir_ol  , '/obs_epoch_001.nc'               ));
    end
end

% Figure the links and time variables in the netcdf file
ncid  = netcdf.open(nc(e).state_mean_pr, 'NC_NOWRITE');
ncvar = netcdf.inqDim(ncid, 0);
if strcmp(ncvar, 'links')
    iL = 0; iT = 1;
else
    iT = 0; iL = 1;
end

% Retrieve the dimensions from the netcdf file
Nt_tmp = zeros(1, num_exps);
for e = 1:num_exps
    ncid            = netcdf.open(nc(e).state_mean_pr, 'NC_NOWRITE');
    [~, Nt_tmp(e)]  = netcdf.inqDim(ncid, iT); % # of assim cycles
    netcdf.close(ncid);
end
if length(unique(Nt_tmp)) > 1
    error([ 'Numer of DA cycles in the ' num2str(num_exps) ' experiments is not the same! Exiting ...' ])
end
Nt = Nt_tmp(1);

ncid    = netcdf.open(nc(1).state_mean_pr, 'NC_NOWRITE');
[~, Nl] = netcdf.inqDim(ncid, iL); % # of links in the domain             
netcdf.close(ncid);

exp = struct;
for e = 1:num_exps
    ncid        = netcdf.open(nc(e).obs_diag, 'NC_NOWRITE');
    [~, rbins]  = netcdf.inqDim(ncid, 14); % # of bins in the rank histogram (i.e., ens_size+1) 

    exp(e).ens_size = rbins-1; % size of the ensemble
    netcdf.close(ncid);
end 

% ensemble size of the open loop
ncid        = netcdf.open(ol.obs_diag, 'NC_NOWRITE');
[~, rbins]  = netcdf.inqDim(ncid, 14); % # of bins in the rank histogram (i.e., ens_size+1) 
ol.ens_size = rbins-1; % size of the ensemble
netcdf.close(ncid);

% separate exp name from path
exp_name = string(missing);
for e = 1:num_exps
    sp_names    = strsplit(dir_exps{e}, '/');
    exp_name(e) = sp_names{end};
end


%% TIME HANDLING:
Time        = double(ncread(nc(1).state_mean_pr, 'time')); 
unit_time   = ncreadatt(nc(1).state_mean_pr,'time','units');

if Time(1) == 0
    Time = Time/24;
end

origin      = datenum(unit_time(12:22));
current     = Time + origin; 
goodtime    = datestr(current, 'mmm dd'); 
longtime    = datestr(current, 'mmm dd, yyyy HH:MM pm');

time_label  = ceil( [1, Nt/4, Nt/2, 3*Nt/4, Nt] );
xticks      = current(time_label);
xtickslabel = goodtime(time_label, :);

% obs_diag time
od_time     = ncread(nc(e).obs_diag, 'time') + origin;
od_time_b1  = find(od_time == current(1));
od_time_b2  = find(od_time == current(Nt));
diag_range  = od_time_b1:od_time_b2; 

%% PROCESSING DATA:
gauge_id = strtrim(ncread(nc(1).routelink, 'gages')');

obserr_L = 1.e8; 

% All available gauges in the domain: 
k = 0;
for i = 1:Nl
    ob_id = gauge_id(i, :);
    if sum( isspace(ob_id) ) < 10 
        k = k+1;
        gauges.avail.OID(k) = str2double(ob_id); % Available gauges IDs
        gauges.avail.IND(k) = i;                 % Available gauges indices
    end 
end
gauges.avail.num  = length(gauges.avail.OID);    % # of all available gauges

gauges.yaml.names = strtrim(ncread(nc(1).obs_diag, 'ObservationTypes')');

for l = 1:size(gauges.yaml.names, 1)
    gauges.yaml.chari(l, :) = str2double(gauges.yaml.names(l, 13:end));
end

if obs == 0
    % display all gauges (wanted in yaml)
    gauges.want.names = gauges.yaml.names;
    gauges.want.num   = size(gauges.want.names, 1);

    gauges.want.IND = zeros(1, gauges.want.num);
    for l = 1:gauges.want.num
        gauges.want.IND(l) = gauges.yaml.chari(l, :);
        gauges.want.OID(l) = gauges.avail.OID(gauges.avail.IND == gauges.want.IND(l));
    end
    
else
    % user-specified gauges
    gauges.want.num = length(obs);
    gauges.want.OID = obs;
    
    for l = 1:gauges.want.num
        gauges.want.IND(l)      = gauges.avail.IND(gauges.avail.OID == gauges.want.OID(l));
        gauges.want.names(l, :) = gauges.yaml.names(gauges.yaml.chari == gauges.want.IND(l), :);
    end
    
end

% Reading:
for e = 1:num_exps
    exp(e).ensemble = double(ncread(nc(e).obs_epoc, 'observations'));
    exp(e).obs_ind  = -1 * double(ncread(nc(e).obs_epoc, 'obs_type'));
    exp(e).O_time   = double(ncread(nc(e).obs_epoc, 'time')) + origin;
    
    if plot_ol
        ol.ensemble = double(ncread(ol.obs_epoc, 'observations')); 
        ol.obs_ind  = -1 * double(ncread(ol.obs_epoc, 'obs_type'));
    end
    
    % State and spread files
    exp(e).pr.state.x1  = ncread(nc(e).state_mean_pr, 'qlink1'); 
    exp(e).pr.spread.x1 = ncread(nc(e).state_sd_pr,   'qlink1');
    exp(e).po.state.x1  = ncread(nc(e).state_mean_po, 'qlink1');
    exp(e).po.spread.x1 = ncread(nc(e).state_sd_po,   'qlink1');
    
    [~, numvars] = netcdf.inq(netcdf.open(nc(e).state_mean_pr, 'NC_NOWRITE'));
    if numvars > 2
        bucket = true;
    else
        bucket = false;
    end
    
    if bucket
        exp(e).pr.state.x2  = ncread(nc(e).state_mean_pr, 'z_gwsubbas') * 1000; % m -> mm 
        exp(e).pr.spread.x2 = ncread(nc(e).state_sd_pr,   'z_gwsubbas') * 1000; % m -> mm 
        exp(e).po.state.x2  = ncread(nc(e).state_mean_po, 'z_gwsubbas') * 1000; % m -> mm 
        exp(e).po.spread.x2 = ncread(nc(e).state_sd_po,   'z_gwsubbas') * 1000; % m -> mm 
    end
    
    % inflation files
    if inf_flav_n(e, 1) > 0 || inf_flav_n(e, 2) > 0
        if inf_flav_n(e, 1) > 0 && inf_flav_n(e, 2) > 0 
            exp(e).pr.infm.x1 = ncread(nc(e).pr_inflate_mean, 'qlink1');
            exp(e).pr.infs.x1 = ncread(nc(e).pr_inflate_sd  , 'qlink1');
            exp(e).po.infm.x1 = ncread(nc(e).po_inflate_mean, 'qlink1');
            exp(e).po.infs.x1 = ncread(nc(e).po_inflate_sd  , 'qlink1');

            if bucket
                exp(e).pr.infm.x2 = ncread(nc(e).pr_inflate_mean, 'z_gwsubbas');
                exp(e).pr.infs.x2 = ncread(nc(e).pr_inflate_sd  , 'z_gwsubbas');
                exp(e).po.infm.x2 = ncread(nc(e).po_inflate_mean, 'z_gwsubbas');
                exp(e).po.infs.x2 = ncread(nc(e).po_inflate_sd  , 'z_gwsubbas');
            end
            
        elseif inf_flav_n(e, 1) > 0 
            exp(e).pr.infm.x1 = ncread(nc(e).pr_inflate_mean, 'qlink1');
            exp(e).pr.infs.x1 = ncread(nc(e).pr_inflate_sd  , 'qlink1');

            if bucket
                exp(e).pr.infm.x2 = ncread(nc(e).pr_inflate_mean, 'z_gwsubbas');
                exp(e).pr.infs.x2 = ncread(nc(e).pr_inflate_sd  , 'z_gwsubbas');
            end
            
        elseif inf_flav_n(e, 2) > 0 
            exp(e).po.infm.x1 = ncread(nc(e).po_inflate_mean, 'qlink1');
            exp(e).po.infs.x1 = ncread(nc(e).po_inflate_sd  , 'qlink1');

            if bucket
                exp(e).po.infm.x2 = ncread(nc(e).po_inflate_mean, 'z_gwsubbas');
                exp(e).po.infs.x2 = ncread(nc(e).po_inflate_sd  , 'z_gwsubbas');
            end
        end

        % Hybrid files
        if hyb_flav_n(e) > 0 
            exp(e).pr.hybm.x1 = ncread(nc(e).hybrid_mean, 'qlink1');
            exp(e).pr.hybs.x1 = ncread(nc(e).hybrid_sd  , 'qlink1');

            if bucket
                exp(e).pr.hybm.x2 = ncread(nc(e).hybrid_mean, 'z_gwsubbas');
                exp(e).pr.hybs.x2 = ncread(nc(e).hybrid_sd  , 'z_gwsubbas');
            end
        end
    end
end

openloop    = cell( 4, gauges.want.num          );
forecast    = cell(12, gauges.want.num, num_exps);
observation = cell( 9, gauges.want.num, num_exps);
analysis    = cell( 9, gauges.want.num, num_exps);

flood       = zeros(gauges.want.num, num_exps);

for i = 1:gauges.want.num

    k = gauges.want.IND(i);
    
    for e = 1:num_exps
    
        find_obs = k == exp(e).obs_ind;

        fprintf('exp: %2d, obs no: %3d,  dart-index: %6d,  USGS gauge ID: %10d\n', e, i, k, gauges.want.OID(i))

        tmp.obs_val    = exp(e).ensemble(1,          find_obs);
        tmp.ens_mean_f = exp(e).ensemble(2,          find_obs);
        tmp.ens_sd_f   = exp(e).ensemble(4,          find_obs);
        tmp.ensemble_f = exp(e).ensemble(6:2:end-1,  find_obs);
        tmp.obs_var    = exp(e).ensemble(end,        find_obs); 
        tmp.ens_mean_a = exp(e).ensemble(3,          find_obs);
        tmp.ens_sd_a   = exp(e).ensemble(5,          find_obs);
        tmp.ensemble_a = exp(e).ensemble(7:2:end,    find_obs);
        
        if plot_ol
            % Just in case, open loop has more gauges
            ol_obs = k == ol.obs_ind;

            tmp.ens_mean_ol = ol.ensemble(2,         ol_obs); 
            tmp.ens_sd_ol   = ol.ensemble(4,         ol_obs); 
            tmp.ensemble_ol = ol.ensemble(6:2:end-1, ol_obs);
        end  

        ens_time   = zeros(1, Nt+1);
        Found_time = exp(e).O_time(find_obs); 
        
        % initial time management:
        % some obs are off by 5 mins
        z_offset = current(1) - Found_time(1);
        if z_offset ~= 0
            Found_time = Found_time + z_offset;
        end
    
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
        ensemble_f = NaN(exp(e).ens_size, Nt);
        obs_var    = zeros(1, Nt);
        ens_mean_a = zeros(1, Nt);
        ens_sd_a   = zeros(1, Nt);
        ensemble_a = NaN(exp(e).ens_size, Nt);
        
        if plot_ol
            ens_mean_ol = zeros(1, Nt); 
            ens_sd_ol   = zeros(1, Nt); 
            ensemble_ol = NaN(ol.ens_size, Nt);
        end
    
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

                if plot_ol
                    ens_mean_ol(j)    = mean(tmp.ens_mean_ol(:, ens_time(j):ens_time(j+1)), 2); 
                    ens_sd_ol(j)      = mean(tmp.ens_sd_ol(:, ens_time(j):ens_time(j+1)), 2); 
                    ensemble_ol(:, j) = mean(tmp.ensemble_ol(:, ens_time(j):ens_time(j+1)), 2);
                end

            elseif ~isnan(ens_time(j))
                obs_val(j)       = tmp.obs_val   (:, ens_time(j));
                ens_mean_f(j)    = tmp.ens_mean_f(:, ens_time(j));
                ens_sd_f(j)      = tmp.ens_sd_f  (:, ens_time(j));
                ensemble_f(:, j) = tmp.ensemble_f(:, ens_time(j));
                obs_var(j)       = tmp.obs_var   (:, ens_time(j));
                ens_mean_a(j)    = tmp.ens_mean_a(:, ens_time(j));
                ens_sd_a(j)      = tmp.ens_sd_a  (:, ens_time(j));
                ensemble_a(:, j) = tmp.ensemble_a(:, ens_time(j));
                
                if plot_ol
                    ens_mean_ol(j)    = tmp.ens_mean_ol(:, ens_time(j));
                    ens_sd_ol(j)      = tmp.ens_sd_ol(:, ens_time(j));
                    ensemble_ol(:, j) = tmp.ensemble_ol(:, ens_time(j));
                end
                    

            elseif ~isnan(ens_time(j+1))
                obs_val(j)       = tmp.obs_val   (:, ens_time(j+1));
                ens_mean_f(j)    = tmp.ens_mean_f(:, ens_time(j+1));
                ens_sd_f(j)      = tmp.ens_sd_f  (:, ens_time(j+1));
                ensemble_f(:, j) = tmp.ensemble_f(:, ens_time(j+1));
                obs_var(j)       = tmp.obs_var   (:, ens_time(j+1));
                ens_mean_a(j)    = tmp.ens_mean_a(:, ens_time(j+1));
                ens_sd_a(j)      = tmp.ens_sd_a  (:, ens_time(j+1));
                ensemble_a(:, j) = tmp.ensemble_a(:, ens_time(j+1));
                
                if plot_ol
                    ens_mean_ol(j)    = tmp.ens_mean_ol(:, ens_time(j+1));
                    ens_sd_ol(j)      = tmp.ens_sd_ol(:, ens_time(j+1));
                    ensemble_ol(:, j) = tmp.ensemble_ol(:, ens_time(j+1));
                end

            else
                obs_val(j)       = NaN;
                ens_mean_f(j)    = NaN;
                ens_sd_f(j)      = NaN;
                ensemble_f(:, j) = NaN;
                obs_var(j)       = NaN;
                ens_mean_a(j)    = NaN;
                ens_sd_a(j)      = NaN;
                ensemble_a(:, j) = NaN;
                
                if plot_ol
                    ens_mean_ol(j)    = NaN; 
                    ens_sd_ol(j)      = NaN; 
                    ensemble_ol(:, j) = NaN;
                end

            end
        end
        clear tmp

        flood(i, e) = find(obs_val == max(obs_val, [], 'omitnan'), 1);
        
        varname_f = strcat(gauges.want.names(i, :), '_guess');
        varname_a = strcat(gauges.want.names(i, :), '_analy');
    
        % Manage open loop copies
        if plot_ol, rmse_ol = abs(ens_mean_ol - obs_val); end
    
        % Manage the forecast copies
        tmp_f        = squeeze(double(ncread(nc(e).obs_diag, varname_f)));
        tmp_f        = tmp_f(:, diag_range);
        rmse_f       = abs(ens_mean_f - obs_val);
        bias_f       = obs_val - ens_mean_f;
        totspread_f  = sqrt(ens_sd_f.^2 + obs_var.^2);
        rank_hist_f  = squeeze(double(ncread(nc(e).obs_diag, strcat(varname_f, '_RankHist'))));
        rank_hist_f  = rank_hist_f(:, diag_range);
        
        if inf_flav_n(e, 1) > 0
            pr_inflate_mean = exp(e).pr.infm.x1(k, :);
            pr_inflate_sd   = exp(e).pr.infs.x1(k, :); 
        elseif inf_flav_n(e, 1) == 0
            pr_inflate_mean = nan;
            pr_inflate_sd   = nan;
        end     
        
        % Save the un-inflated prior
        ensemble_f_def = ensemble_f;
        if inf_flav_n(e, 1) > 0
            lambda = pr_inflate_mean;
            for j = 1:Nt %deflate the ensemble
                ensemble_f_def(:, j) = 1/lambda(j) * (ensemble_f(:, j) - ens_mean_f(j)) + ens_mean_f(j);
            end
        end

        if hyb_flav_n(e) > 0
            hybrid_mean = exp(e).pr.hybm.x1(k, :);
            hybrid_sd   = exp(e).pr.hybs.x1(k, :); 
        else
            hybrid_mean = nan;
            hybrid_sd   = nan;
        end
    
        % Manage the observation copies
        obs_poss    = tmp_f(1, :);
        obs_used    = tmp_f(2, :);
        obs_w_qc7   = tmp_f(22, :);
    
        obs_assimilated = obs_used ~= 0;
        obs_rejected    = obs_used == 0;
        obs_val_assim   = obs_val;
        obs_val_reject  = obs_val;

        obs_val_assim(obs_rejected)     = NaN;
        obs_val_reject(obs_assimilated) = NaN;
    
        % Manage the analysis copies
        rmse_a      = abs(ens_mean_a - obs_val);
        bias_a      = obs_val - ens_mean_a;
        totspread_a = sqrt(ens_sd_a.^2 + obs_var.^2);
        
        if inf_flav_n(e, 2) > 0
            po_inflate_mean = exp(e).po.infm.x1(k, :);
            po_inflate_sd   = exp(e).po.infs.x1(k, :);
            
        elseif inf_flav_n(e, 2) == 0
            po_inflate_mean = nan;
            po_inflate_sd   = nan;
        end 
        
        % Save the un-inflated prior
        ensemble_a_def = ensemble_a;
        if inf_flav_n(e, 2) > 0
            lambda = po_inflate_mean;
            for j = 1:Nt %deflate the ensemble
                ensemble_a_def(:, j) = 1/lambda(j) * (ensemble_a(:, j) - ens_mean_a(j)) + ens_mean_a(j);
            end
        end 
    
    
        % Construct the (open loop), forecast, observation and analysis cells.
        if plot_ol
            openloop(:, i)       = { ens_mean_ol    ; ...
                                     rmse_ol        ; ...
                                     ens_sd_ol      ; ...
                                     ensemble_ol    ; ...
                                   };
        end
        
        forecast(:, i, e)    = { varname_f       ; ...
                                 rmse_f          ; ...
                                 bias_f          ; ... 
                                 ens_sd_f        ; ...
                                 totspread_f     ; ...
                                 ens_mean_f      ; ...
                                 rank_hist_f     ; ...
                                 ensemble_f_def  ; ...
                                 pr_inflate_mean ; ...
                                 pr_inflate_sd   ; ...
                                 hybrid_mean     ; ...
                                 hybrid_sd       ; ...
                               };
                 
        observation(:, i, e) = { k                  ; ...
                                 gauges.want.OID(i) ; ...
                                 obs_poss           ; ...
                                 obs_used           ; ...
                                 obs_w_qc7          ; ...
                                 obs_val            ; ...
                                 obs_val_assim      ; ...
                                 obs_val_reject     ; ...
                                 obs_var            ; ...
                               };  
                          
                 
        analysis(:, i, e)    = { varname_a       ; ...
                                 rmse_a          ; ...
                                 bias_a          ; ... 
                                 ens_sd_a        ; ...
                                 totspread_a     ; ...
                                 ens_mean_a      ; ...
                                 ensemble_a_def  ; ...
                                 po_inflate_mean ; ...
                                 po_inflate_sd   ; ...
                               };
    end
end


%% DISPLAY RESULTS: %
% ***************** %

if disp_res

if isfile(fig_to_pdf), delete(fig_to_pdf); end

%% TIME SERIES EVOLUTION:

for o = 1:gauges.want.num
    
    if num_exps > 1
        
        if num_exps == 2
            figure('uni','pi','pos',[50, 600, 1600, 450]); 
        elseif num_exps <= 4
            figure('uni','pi','pos',[50, 600, 1600, 900]);
        else
            figure('uni','pi','pos',[50, 600, 1600, 1000]);
        end
            
        
        for e = 1:num_exps
            
            if num_exps == 2

                start_x = [.05, .55];
                fig_wid = .40;
                
                fig_ht1 = .50;
                fig_ht2 = .26;
                
                sep     = .05;
                bot_y   = .08;
                top_y   = bot_y + sep + fig_ht2;

                subplot('Position', [start_x(e), top_y, fig_wid, fig_ht1]);
                
                en_f = plot(current, forecast{8, o, e}   , '-', 'Color', gY); hold on
                en_a = plot(current, analysis{7, o, e}   , '-', 'Color', lB); grid on 
                
                ob_a = plot(current, observation{7, o, e}, '*', 'Color', gR); 
                ob_r = plot(current, observation{8, o, e}, '*', 'Color', rD); 
                
                if plot_ol, op   = plot(current, openloop{1, o}      , '-', 'Color', oR, 'LineWidth', 3); end
                mf   = plot(current, forecast{6, o, e}   , '-', 'Color', bK, 'LineWidth', 3); 
                ma   = plot(current, analysis{6, o, e}   , '-', 'Color', bL, 'LineWidth', 3); 

                if plot_ol, so   = plot(current, openloop{3, o}      , '--', 'Color', oR, 'LineWidth', 1); end
                sf   = plot(current, forecast{4, o, e}   , '--', 'Color', bK, 'LineWidth', 1); 
                sa   = plot(current, analysis{4, o, e}   , '--', 'Color', bL, 'LineWidth', 1); 
                
                limsy = get(gca, 'YLim');
                
                strY = 'Sream flow (cms)'; 

                set(gca, 'FontSize', 16, 'XLim', [xticks(1), xticks(end)], 'XTick', xticks, 'XTickLabel', {}, 'Ylim', [0 limsy(2)])
                ylabel(strY, 'FontSize', 18)

                if mean(observation{9, o, e}, 'omitnan') < obserr_L
                    obs_status = 'Used Obs (Assimilated)';
                else
                    obs_status = 'Used Obs (Evaluated only)';
                end

                if plot_ol
                    L = legend([ob_a, ob_r, en_f(1), en_a(1), op, mf , ma, so, sf, sa], ...
                                obs_status          , ...
                                sprintf('Rejected Obs: %.2f%%'       , 100-observation{4, o, e}/observation{3, o, e}*100), ...
                                'Prior Members', 'Posterior Members' , ...
                                sprintf('Open Loop, RMSE: %.2f'      , mean(openloop{2, o}, 'omitnan')), ...
                                sprintf('Prior Mean, RMSE: %.2f'     , mean(forecast{2, o, e}, 'omitnan')), ...
                                sprintf('Posterior Mean, RMSE: %.2f' , mean(analysis{2, o, e}, 'omitnan')), ...
                                sprintf('Open loop Spread, avg: %.2f', mean(openloop{3, o}, 'omitnan')), ... 
                                sprintf('Prior Spread, avg: %.2f'    , mean(forecast{4, o, e}, 'omitnan')), ...
                                sprintf('Posterior Spread, avg: %.2f', mean(analysis{4, o, e}, 'omitnan')), ...
                                'Location', 'NorthEast');
                else
                            L = legend([ob_a, ob_r, mf , ma], obs_status          , ...
                                sprintf('Rejected Obs: %.2f%%'       , 100-observation{4, o, e}/observation{3, o, e}*100), ...
                                sprintf('Prior Mean, RMSE: %.2f'     , mean(forecast{2, o, e}, 'omitnan')), ...
                                sprintf('Posterior Mean, RMSE: %.2f' , mean(analysis{2, o, e}, 'omitnan')), ...
                                'Location', 'NorthEast');
                end
                        
                set(L, 'Interpreter', 'none', 'FontSize', 12, 'color', 'none') 
                title(L, exp_name(e), 'FontSize', 14)

                str2 = [ 'Hydrograph: ', gauges.want.names(o, :), ', Gauge ID: ', num2str(gauges.want.OID(o)) ];
                title(str2, 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none')

                subplot('Position', [start_x(e), bot_y, fig_wid, fig_ht2]);

                if sum(inf_flav_n) > 1
                    i_pr_m = plot(current, forecast{9, o, e}, '-', 'Color', bK, 'LineWidth', 2); hold on
                    i_po_m = plot(current, analysis{8, o, e}, '-', 'Color', bL, 'LineWidth', 2); grid on
                    
                    plot(current, ones(1, Nt), '--', 'Color', gY); ax1 = gca;
                    
                    limsy = get(gca, 'YLim');
                    
                    set(gca, 'FontSize', 16, 'XLim', [xticks(1), xticks(end)], 'XTick', xticks, 'XTickLabel', xtickslabel, 'Ylim', [0.9 limsy(2)])
                    ylabel('Inflation', 'FontSize', 18)
                    
                    if hyb_flav_n(e) > 0
                        yyaxis right
                        hyb_m = plot(current, forecast{11, o, e}, '-', 'Color', gR, 'LineWidth', 2); 
                        set(gca, 'YColor', bK, 'YLim', [-0.1, 1.1], 'YTick', [0, 0.5, 1])
                        ylabel('Hyb. Weight', 'FontSize', 18)

                        L = legend(ax1, [i_pr_m, i_po_m, hyb_m], ...
                            sprintf('Prior Inflation Mean, avg: %.2f'      , mean(forecast{9,  o, e}, 'omitnan')), ...
                            sprintf('Posterior Inflation Mean, avg: %.2f'  , mean(analysis{8,  o, e}, 'omitnan')), ...
                            sprintf('Hybrid Weight Mean, avg: %.2f'        , mean(forecast{11,  o, e}, 'omitnan')), ...
                            'Location', 'NorthEast');
                    else
                        L = legend(ax1, [i_pr_m, i_po_m], ...
                            sprintf('Prior Inflation Mean, avg: %.2f'      , mean(forecast{9,  o, e}, 'omitnan')), ...
                            sprintf('Posterior Inflation Mean, avg: %.2f'  , mean(analysis{8,  o, e}, 'omitnan')), ...
                            'Location', 'NorthEast');
                    end
            
                    set(L, 'Interpreter', 'none', 'Box', 'off', 'FontSize', 12, 'color', 'none')
                else
                    i_pr_m = plot(current, forecast{9, o, e}, '-', 'Color', bK, 'LineWidth', 2); hold on
                    
                    plot(current, ones(1, Nt), '--', 'Color', gY); grid on
                    
                    limsy = get(gca, 'YLim');
                    
                    set(gca, 'FontSize', 16, 'XLim', [xticks(1), xticks(end)], 'XTick', xticks, 'XTickLabel', xtickslabel, 'Ylim', [0.9 limsy(2)])
                    ylabel('Inflation', 'FontSize', 18)

                    if hyb_flav_n(e) > 0
                        yyaxis right
                        hyb_m = plot(current, forecast{11, o, e}, '-', 'Color', gR, 'LineWidth', 2); 
                        set(gca, 'YColor', bK, 'YLim', [-0.1, 1.1], 'YTick', [0, 0.5, 1])
                        ylabel('Hyb. Weight', 'FontSize', 18)

                        L = legend(ax1, [i_pr_m, hyb_m], ...
                            sprintf('Prior Inflation Mean, avg: %.2f'      , mean(forecast{9,  o, e}, 'omitnan')), ...
                            sprintf('Hybrid Weight Mean, avg: %.2f'        , mean(forecast{11,  o, e}, 'omitnan')), ...
                            'Location', 'NorthEast');
                    else
                        L = legend(sprintf('Prior Inflation Mean, avg: %.2f', mean(forecast{9,  o, e}, 'omitnan')), ...
                            'Location', 'NorthEast');
                    end
            
                    set(L, 'Interpreter', 'none', 'Box', 'off', 'FontSize', 12, 'color', 'none')
                end
               
            else
                
                % num_exps > 2 :: no place for inflation!
                rows = ceil(num_exps/2);
                subplot(rows, 2, e)
                
                en_f = plot(current, forecast{8, o, e}   , '-', 'Color', gY); hold on
                en_a = plot(current, analysis{7, o, e}   , '-', 'Color', lB); grid on 
                
                ob_a = plot(current, observation{7, o, e}, '*', 'Color', gR); 
                ob_r = plot(current, observation{8, o, e}, '*', 'Color', rD);
                
                if plot_ol, op   = plot(current, openloop{1, o}      , '-', 'Color', oR, 'LineWidth', 3); end
                mf   = plot(current, forecast{6, o, e}   , '-', 'Color', bK, 'LineWidth', 3); 
                ma   = plot(current, analysis{6, o, e}   , '-', 'Color', bL, 'LineWidth', 3); 

                if plot_ol, so   = plot(current, openloop{3, o}      , '--', 'Color', oR, 'LineWidth', 1); end
                sf   = plot(current, forecast{4, o, e}   , '--', 'Color', bK, 'LineWidth', 1); 
                sa   = plot(current, analysis{4, o, e}   , '--', 'Color', bL, 'LineWidth', 1); 

                limsy = get(gca, 'YLim');

                set(gca, 'FontSize', 16, 'XLim', [xticks(1), xticks(end)], 'XTick', xticks, 'XTickLabel', xtickslabel, 'Ylim', [0 limsy(2)])
                ylabel(['Gauge: ', num2str(observation{2, o, e})], 'FontSize', 18)

                if mean(observation{9, o, e}, 'omitnan') < obserr_L
                    obs_status = 'Used Obs (Assimilated)';
                else
                    obs_status = 'Used Obs (Evaluated only)';
                end

                if plot_ol
                    L = legend([ob_a, ob_r, en_f(1), en_a(1), op, mf , ...
                                ma, so, sf, sa], obs_status          , ...
                                sprintf('Rejected Obs: %.2f%%'       , 100-observation{4, o, e}/observation{3, o, e}*100), ...
                                'Prior Members', 'Posterior Members' , ...
                                sprintf('Open Loop, RMSE: %.2f'      , mean(openloop{2, o}, 'omitnan')), ...
                                sprintf('Prior Mean, RMSE: %.2f'     , mean(forecast{2, o, e}, 'omitnan')), ...
                                sprintf('Posterior Mean, RMSE: %.2f' , mean(analysis{2, o, e}, 'omitnan')), ...
                                sprintf('Open loop Spread, avg: %.2f', mean(openloop{3, o}, 'omitnan')), ... 
                                sprintf('Prior Spread, avg: %.2f'    , mean(forecast{4, o, e}, 'omitnan')), ...
                                sprintf('Posterior Spread, avg: %.2f', mean(analysis{4, o, e}, 'omitnan')), ...
                                'Location', 'NorthEast');
                else
                            L = legend([ob_a, ob_r, en_f(1), en_a(1), mf , ...
                                ma, sf, sa], obs_status          , ...
                                sprintf('Rejected Obs: %.2f%%'       , 100-observation{4, o, e}/observation{3, o, e}*100), ...
                                'Prior Members', 'Posterior Members' , ...
                                sprintf('Prior Mean, RMSE: %.2f'     , mean(forecast{2, o, e}, 'omitnan')), ...
                                sprintf('Posterior Mean, RMSE: %.2f' , mean(analysis{2, o, e}, 'omitnan')), ...
                                sprintf('Prior Spread, avg: %.2f'    , mean(forecast{4, o, e}, 'omitnan')), ...
                                sprintf('Posterior Spread, avg: %.2f', mean(analysis{4, o, e}, 'omitnan')), ...
                                'Location', 'NorthEast');
                end                        
                set(L, 'Interpreter', 'none', 'FontSize', 10, 'color', 'none') 

                title('Exp: ' + exp_name(e), 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none')
                
            end
                
        end
        if fig_to_pdf
            exportgraphics(gcf, fig_to_pdf, 'Append', true, 'ContentType', 'vector')
            close
        end
        
    else
    
        % Only 1 experiment
        figure('uni','pi','pos',[200, 600, 800, 500]) 
        
        start_x = .10;
        fig_wid = .80;
        
        fig_ht1 = .50;
        fig_ht2 = .26;
        
        sep     = .05;
        bot_y   = .08;
        top_y   = bot_y + sep + fig_ht2;
       
        
        subplot('Position', [start_x, top_y, fig_wid, fig_ht1]);
        
        en_f = plot(current, forecast{8, o, e}   , '-', 'Color', gY); hold on
        en_a = plot(current, analysis{7, o, e}   , '-', 'Color', lB); grid on
        
        ob_a = plot(current, observation{7, o, e}, '*', 'Color', gR);  
        ob_r = plot(current, observation{8, o, e}, '*', 'Color', rD); 
        
        if plot_ol, op   = plot(current, openloop{1, o}      , '-', 'Color', oR, 'LineWidth', 3); end
        mf   = plot(current, forecast{6, o, e}   , '-', 'Color', bK, 'LineWidth', 3); 
        ma   = plot(current, analysis{6, o, e}   , '-', 'Color', bL, 'LineWidth', 3); 

        if plot_ol, so   = plot(current, openloop{3, o}      , '--', 'Color', oR, 'LineWidth', 1); end
        sf   = plot(current, forecast{4, o, e}   , '--', 'Color', bK, 'LineWidth', 1); 
        sa   = plot(current, analysis{4, o, e}   , '--', 'Color', bL, 'LineWidth', 1); 
        
        limsy = get(gca, 'YLim');

        set(gca, 'FontSize', 16, 'XLim', [xticks(1), xticks(end)], 'XTick', xticks, 'XTickLabel', {}, 'Ylim', [0 limsy(2)])
        ylabel('Stream flow (cms)', 'FontSize', 18)

        if mean(observation{9, o, e}, 'omitnan') < obserr_L
            obs_status = 'Used Obs (Assimilated)';
        else
            obs_status = 'Used Obs (Evaluated only)';
        end

        if plot_ol 
            L = legend([ob_a, ob_r, en_f(1), en_a(1), op, mf, ma, ...
                    so, sf, sa], obs_status              , ...
                    sprintf('Rejected Obs: %.2f%%'       , 100-observation{4, o, e}/observation{3, o, e}*100), ...
                    'Prior Members', 'Posterior Members' , ...
                    sprintf('Open Loop, RMSE: %.2f'      , mean(openloop{2,  o   }, 'omitnan')), ...
                    sprintf('Prior Mean, RMSE: %.2f'     , mean(forecast{2,  o, e}, 'omitnan')), ...
                    sprintf('Posterior Mean, RMSE: %.2f' , mean(analysis{2,  o, e}, 'omitnan')), ...
                    sprintf('Open loop Spread, avg: %.2f', mean(openloop{3,  o   }, 'omitnan')), ...
                    sprintf('Prior Spread, avg: %.2f'    , mean(forecast{4,  o, e}, 'omitnan')), ...
                    sprintf('Posterior Spread, avg: %.2f', mean(analysis{4,  o, e}, 'omitnan')), ...
                    'Location', 'NorthEast');
        else
            L = legend([ob_a, ob_r, en_f(1), en_a(1), mf, ma, ...
                    sf, sa], obs_status              , ...
                    sprintf('Rejected Obs: %.2f%%'       , 100-observation{4, o, e}/observation{3, o, e}*100), ...
                    'Prior Members', 'Posterior Members' , ...
                    sprintf('Prior Mean, RMSE: %.2f'     , mean(forecast{2,  o, e}, 'omitnan')), ...
                    sprintf('Posterior Mean, RMSE: %.2f' , mean(analysis{2,  o, e}, 'omitnan')), ...
                    sprintf('Prior Spread, avg: %.2f'    , mean(forecast{4,  o, e}, 'omitnan')), ...
                    sprintf('Posterior Spread, avg: %.2f', mean(analysis{4,  o, e}, 'omitnan')), ...
                    'Location', 'NorthEast');
        end
        set(L, 'Interpreter', 'none', 'Box', 'off', 'FontSize', 12, 'color', 'none')

        str1 = [ gauges.want.names(o, :), ' | Gauge ID: ', num2str(gauges.want.OID(o)) ];
        str2 = 'Hydrograph: Obs, Prior/Posterior Ensemble, Mean, Spread & Inflation';
        title({str1, str2}, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none')
        
        subplot('Position', [start_x, bot_y, fig_wid, fig_ht2]);
        
        if sum(inf_flav_n) > 1
            i_pr_m = plot(current, forecast{9, o, e}, '-', 'Color', bK, 'LineWidth', 2); hold on
            i_po_m = plot(current, analysis{8, o, e}, '-', 'Color', bL, 'LineWidth', 2); grid on
            
            plot(current, ones(1, Nt), '--', 'Color', gY); ax1 = gca;
            
            limsy = get(gca, 'YLim');
            
            set(gca, 'FontSize', 16, 'XLim', [xticks(1), xticks(end)], 'XTick', xticks, 'XTickLabel', xtickslabel, 'Ylim', [0.9 limsy(2)])
            ylabel('Inflation', 'FontSize', 18)

            if hyb_flav_n(e) > 0
                yyaxis right
                hyb_m = plot(current, forecast{11, o, e}, '-', 'Color', gR, 'LineWidth', 2); 
                set(gca, 'YColor', bK, 'YLim', [-0.1, 1.1], 'YTick', [0, 0.5, 1])
                ylabel('Hyb. Weight', 'FontSize', 18)

                L = legend(ax1, [i_pr_m, i_po_m, hyb_m], ...
                    sprintf('Prior Inflation Mean, avg: %.2f'      , mean(forecast{9,  o, e}, 'omitnan')), ...
                    sprintf('Posterior Inflation Mean, avg: %.2f'  , mean(analysis{8,  o, e}, 'omitnan')), ...
                    sprintf('Hybrid Weight Mean, avg: %.2f'        , mean(forecast{11,  o, e}, 'omitnan')), ...
                    'Location', 'NorthEast');
            else
                L = legend(ax1, [i_pr_m, i_po_m], ...
                    sprintf('Prior Inflation Mean, avg: %.2f'      , mean(forecast{9,  o, e}, 'omitnan')), ...
                    sprintf('Posterior Inflation Mean, avg: %.2f'  , mean(analysis{8,  o, e}, 'omitnan')), ...
                    'Location', 'NorthEast');
            end
    
            set(L, 'Interpreter', 'none', 'Box', 'off', 'FontSize', 12, 'color', 'none')
        else
            i_pr_m = plot(current, forecast{9, o, e}, '-', 'Color', bK, 'LineWidth', 2); hold on
            
            plot(current, ones(1, Nt), '--', 'Color', gY); grid on
            
            limsy = get(gca, 'YLim');
            
            set(gca, 'FontSize', 16, 'XLim', [xticks(1), xticks(end)], 'XTick', xticks, 'XTickLabel', xtickslabel, 'Ylim', [0.9 limsy(2)])
            ylabel('Inflation', 'FontSize', 18)

            if hyb_flav_n(e) > 0
                yyaxis right
                hyb_m = plot(current, forecast{11, o, e}, '-', 'Color', gR, 'LineWidth', 2); 
                set(gca, 'YColor', bK, 'YLim', [-0.1, 1.1], 'YTick', [0, 0.5, 1])
                ylabel('Hyb. Weight', 'FontSize', 18)

                L = legend(ax1, [i_pr_m, hyb_m], ...
                    sprintf('Prior Inflation Mean, avg: %.2f'      , mean(forecast{9,  o, e}, 'omitnan')), ...
                    sprintf('Hybrid Weight Mean, avg: %.2f'        , mean(forecast{11,  o, e}, 'omitnan')), ...
                    'Location', 'NorthEast');
            else
                L = legend(sprintf('Prior Inflation Mean, avg: %.2f', mean(forecast{9,  o, e}, 'omitnan')), ...
                    'Location', 'NorthEast');
            end
    
            set(L, 'Interpreter', 'none', 'Box', 'off', 'FontSize', 12, 'color', 'none')
        end

        if fig_to_pdf
            %print(gcf, '-dpsc', '-vector', '-append', '-bestfit', pdf_filename)
            exportgraphics(gcf, fig_to_pdf, 'Append', true, 'ContentType', 'vector')
            close
        end
        
        if plot_state && o == gauges.want.num
            
            tiny_flow_s = 10;
            tiny_flow_b = 1;
            
            % display avg. mean and spread
            figure('uni','pi','pos',[200, 600, 1200, 970])
        
            Xm1_f = mean(exp(e).pr.state.x1 , 2); 
            Xs1_f = mean(exp(e).pr.spread.x1, 2); 
            Xm1_a = mean(exp(e).po.state.x1 , 2);  
            Xs1_a = mean(exp(e).po.spread.x1, 2); 

            subplot(221)
            plot_connections(Xm1_f, nc(e).routelink, get(gca, 'position'), 'cms', tiny_flow_s)
            title({'Experiment: ' + exp_name(e),'Stream Flow: Time-Avg. Prior Mean'}, 'FontSize', 14, 'Interpreter', 'none')

            subplot(222)
            plot_connections(Xs1_f, nc(e).routelink, get(gca, 'position'), 'cms', tiny_flow_s)
            title({'Experiment: ' + exp_name(e),'Stream Flow: Time-Avg. Prior Spread'}, 'FontSize', 14, 'Interpreter', 'none')

            subplot(223)
            if bucket
                Xm2_f = mean(exp(e).pr.state.x2 , 2); 
                plot_connections(Xm2_f, nc(e).routelink, get(gca, 'position'), 'mm', tiny_flow_b)
                title({'Experiment: ' + exp_name(e),'Bucket: Time-Avg. Prior Mean'}, 'FontSize', 16, 'Interpreter', 'none')
            else
                plot_connections(Xm1_a, nc(e).routelink, get(gca, 'position'), 'cms', tiny_flow_s)
                title({'Experiment: ' + exp_name(e),'Stream Flow: Time-Avg. Posterior Mean'}, 'FontSize', 16, 'Interpreter', 'none')
            end
                

            subplot(224)
            if bucket
                Xs2_f = mean(exp(e).pr.spread.x2, 2);
                plot_connections(Xs2_f, nc(e).routelink, get(gca, 'position'), 'mm', tiny_flow_b)
                title({'Experiment: ' + exp_name(e),'Bucket: Time-Avg. Prior Spread'}, 'FontSize', 16, 'Interpreter', 'none')
            else
                plot_connections(Xs1_a, nc(e).routelink, get(gca, 'position'), 'cms', tiny_flow_s)
                title({'Experiment: ' + dir_exps(e),'Stream Flow: Time-Avg. Posterior Spread'}, 'FontSize', 16, 'Interpreter', 'none')
            end
            
            % display increment
            figure('uni','pi','pos',[200, 600, 1200, 470])
            
            
            event = ceil(mean(flood(:, e), 1)); 
            Xi1   = exp(e).pr.state.x1(:, event) - exp(e).po.state.x1(:, event);
            Xia   = mean(exp(e).pr.state.x1, 2)  - mean(exp(e).po.state.x1, 2); 
            
            subplot(121)
            plot_connections(Xi1, nc(e).routelink, get(gca, 'position'), 'cms', tiny_flow_s)
            title({'Experiment: ' + exp_name(e),'Stream Flow: DA Increment (Prior-Posterior)', ...
                   ['Event: ' longtime(event, :)]}, 'FontSize', 16, 'Interpreter', 'none')
               
            subplot(122)
            if bucket
                Xi2   = exp(e).pr.state.x2(:, event) - exp(e).po.state.x2(:, event);
                plot_connections(Xi2, nc(e).routelink, get(gca, 'position'), 'mm', tiny_flow_b)
                title({'Experiment: ' + exp_name(e),'Bucket: DA Increment (Prior-Posterior)', ...
                       ['Event: ' longtime(event, :)]}, 'FontSize', 16, 'Interpreter', 'none')
            else
                plot_connections(Xia, nc(e).routelink, get(gca, 'position'), 'cms', tiny_flow_s)
                title({'Experiment: ' + exp_name(e),'Stream Flow: DA Increment (Prior-Posterior)', ...
                   'Time-Average'}, 'FontSize', 16, 'Interpreter', 'none')
            end

        end
 
    end
    
end
      
end 


% % <next few lines under version control, do not edit>
% % $URL: $
% % $Revision: $
% % $Date: $
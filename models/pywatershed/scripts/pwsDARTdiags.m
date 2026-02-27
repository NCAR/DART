function [observation, openloop, forecast, analysis, exp] = ...
    pwsDARTdiags(dir_exps, dir_ol, disp_res, figs_dir, form)

diag_ol = true;

if nargin < 1
    error('No arguments provided!! Please enter the name of experiment directories')
elseif nargin == 1 % no open-loop case
    diag_ol    = false; 
    disp_res   = 0;
elseif nargin <= 3
    disp_res   = 0;
elseif nargin == 4
    form       = '.png';
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
domain_1 = 'd01';
domain_2 = 'd02';

if num_exps > 4
    warning('Too many experiments produce a lot of "cluttered" figures!!')
end

% prior or posterior inflation
inf_flavor = cell(num_exps, 2);
inf_flav_n = zeros(num_exps, 2);

isprior = strcat('/all_output_priorinf_mean_', domain_1, '.nc'); 
ispost  = strcat('/all_output_postinf_mean_', domain_1, '.nc'); 
    
for e = 1:num_exps
    if exist(char(strcat(dir_exps(e), isprior)), 'file') == 2 && ...
        exist(char(strcat(dir_exps(e), ispost)), 'file') == 2
        inf_flavor{e, 1} = '/all_output_priorinf';
        inf_flavor{e, 2} = '/all_output_postinf';
        inf_flav_n(e, 1) = 1;
        inf_flav_n(e, 2) = 1;
        
    elseif exist(char(strcat(dir_exps(e), isprior)), 'file') == 2
        inf_flavor{e, 1} = '/all_output_priorinf'; inf_flavor{e, 2} = '';
        inf_flav_n(e, 1) = 1; inf_flav_n(e, 2) = 0;
        
    elseif exist(char(strcat(dir_exps(e), ispost)), 'file') == 2
        inf_flavor{e, 1} = ''; inf_flavor{e, 2} = '/all_output_postinf';
        inf_flav_n(e, 1) = 0;  inf_flav_n(e, 2) = 1;
    end
end

% Do we have any hybrid weight files?
ishybrid = strcat('/all_output_hybridweight_mean', domain_1, '.nc'); 

hyb_flavor = cell(num_exps);
hyb_flav_n = zeros(num_exps);
for e = 1:num_exps
    if exist(char(strcat(dir_exps(e), ishybrid)), 'file') == 2
        hyb_flavor{e} = '/all_output_hybridweight';
        hyb_flav_n(e) = 1; 
    end
end

nc = struct;
for e = 1:num_exps
    
    diag_dir = dir_exps(e);

    nc(e).routelink     = char(strcat(diag_dir, '/../parameters_dis_seg_app.nc')); % routelink file
    nc(e).obs_diag      = char(strcat(diag_dir, '/obs_diag_output.nc'          )); % output of obs_diag
    nc(e).obs_epoc      = char(strcat(diag_dir, '/obs_epoch_001.nc'            )); % output of obs_seq_to_netcdf
    
    if inf_flav_n(e, 1) > 0 && inf_flav_n(e, 2) > 0
        nc(e).pr_inflate_mean  = char(strcat(diag_dir, inf_flavor{e, 1}, '_mean_', domain_1, '.nc')); % aggregated inf_mean
        nc(e).pr_inflate_sd    = char(strcat(diag_dir, inf_flavor{e, 1}, '_sd_', domain_1, '.nc'  )); % aggregated inf_std 
        nc(e).po_inflate_mean  = char(strcat(diag_dir, inf_flavor{e, 2}, '_mean_', domain_1, '.nc')); % aggregated inf_mean
        nc(e).po_inflate_sd    = char(strcat(diag_dir, inf_flavor{e, 2}, '_sd_', domain_1, '.nc'  )); % aggregated inf_std 
     
    elseif inf_flav_n(e, 1) > 0
        nc(e).pr_inflate_mean  = char(strcat(diag_dir, inf_flavor{e, 1}, '_mean_', domain_1, '.nc')); % aggregated inf_mean
        nc(e).pr_inflate_sd    = char(strcat(diag_dir, inf_flavor{e, 1}, '_sd_', domain_1, '.nc'  )); % aggregated inf_std 
     
    elseif inf_flav_n(e, 2) > 0
        nc(e).po_inflate_mean  = char(strcat(diag_dir, inf_flavor{e, 2}, '_mean_', domain_1, '.nc')); % aggregated inf_mean
        nc(e).po_inflate_sd    = char(strcat(diag_dir, inf_flavor{e, 2}, '_sd_', domain_1, '.nc'  )); % aggregated inf_std 
    end  
 
    if hyb_flav_n(e) > 0
        nc(e).hybrid_mean = char(strcat(diag_dir, hyb_flavor{e}, '_mean_', domain_1, '.nc')); % aggregated hyb_mean
        nc(e).hybrid_sd   = char(strcat(diag_dir, hyb_flavor{e}, '_sd_', domain_1, '.nc'  )); % aggregated hyb_std 
    end 

    % open loop
    if diag_ol 
        ol.obs_diag     = char(strcat(dir_ol  , '/obs_diag_output.nc')); 
        ol.obs_epoc     = char(strcat(dir_ol  , '/obs_epoch_001.nc'  ));
    end
end

% Retrieve the dimensions from the netcdf file

% 1. Number of DA cycles
Nt_tmp = zeros(1, num_exps);
for e = 1:num_exps
    ncid            = netcdf.open(nc(e).obs_diag, 'NC_NOWRITE');
    [~, Nt_tmp(e)]  = netcdf.inqDim(ncid, 0); % # of assim cycles
    netcdf.close(ncid);
end
if length(unique(Nt_tmp)) > 1
    error([ 'Numer of DA cycles in the ' num2str(num_exps) ' experiments is not the same! Exiting ...' ])
end
Nt = Nt_tmp(1);

% 2. Number of links 
ncid    = netcdf.open(nc(1).routelink, 'NC_NOWRITE');
[~, Nl] = netcdf.inqDim(ncid, 0); % # of links in the domain             
netcdf.close(ncid);

% 3. Number of ensemble members
exp = struct;
for e = 1:num_exps
    ncid        = netcdf.open(nc(e).obs_diag, 'NC_NOWRITE');
    [~, rbins]  = netcdf.inqDim(ncid, 14); % # of bins in the rank histogram (i.e., ens_size+1) 

    exp(e).ens_size = rbins-1; % size of the ensemble
    netcdf.close(ncid);
end 

% ensemble size of the open loop
if diag_ol 
    ncid        = netcdf.open(ol.obs_diag, 'NC_NOWRITE');
    [~, rbins]  = netcdf.inqDim(ncid, 14); % # of bins in the rank histogram (i.e., ens_size+1) 
    ol.ens_size = rbins-1; % size of the ensemble
    netcdf.close(ncid);
end

% separate exp name from path
exp_name = string(missing);
for e = 1:num_exps
    sp_names    = strsplit(dir_exps{e}, '/');
    exp_name(e) = sp_names{end};
end

%% TIME HANDLING:
Time        = double(ncread(nc(1).obs_diag, 'time')); 
unit_time   = ncreadatt(nc(1).obs_diag,'time','units');

if Time(1) == 0
    Time = Time/24;
end

origin      = datenum(extractAfter(unit_time,'since '));
current     = Time + origin; 
goodtime    = datestr(current, 'mmm dd'); 
longtime    = datestr(current, 'mmm dd, yyyy HH:MM pm'); %#ok

time_label  = ceil( [1, Nt/4, Nt/2, 3*Nt/4, Nt] );
xticks      = current(time_label);
xtickslabel = goodtime(time_label, :);

% obs_diag time
diag_range  = 1:Nt;

%% PROCESSING DATA:
gauges.avail.OID = str2double(ncread(nc(1).routelink, 'poi_gage_id'));
gauges.avail.IND = double(ncread(nc(1).routelink, 'poi_gage_segment'));

obserr_L = 1.e8; 

% Remove Outside the US gauges
gauges.avail.IND(isnan(gauges.avail.OID)) = [];
gauges.avail.OID(isnan(gauges.avail.OID)) = [];

gauges.avail.num  = length(gauges.avail.OID);    % # of all available gauges

gauges.yaml.names = strtrim(ncread(nc(1).obs_diag, 'ObservationTypes')');

for l = 1:size(gauges.yaml.names, 1)
    gauges.yaml.chari(l, :) = str2double(gauges.yaml.names(l, 13:end));
end

% Display all diagnosed gauges 
X = ncinfo(nc(1).obs_diag);
Z = {X.Attributes(1,:).Name};

diag_gauges = true;
att_count   = length(X.Attributes);
ii = 0;
while(diag_gauges)
    ii = ii + 1;
    try_gauge = str2double(extractAfter(Z{att_count}, 'STREAM_FLOW_'));

    if isnan(try_gauge)
        diag_gauges = false;
    else
        collect_gauge_IDs(ii) = try_gauge; %#ok
    end
    
    att_count = att_count - 1;
end
gauges.want.num = length(collect_gauge_IDs);
state_tag       = 'STREAM_FLOW_';

gauges.want.IND = sort(collect_gauge_IDs);
for l = 1:gauges.want.num
    gauges.want.OID(l) = gauges.avail.OID(gauges.avail.IND == gauges.want.IND(l));
end

% Reading:
for e = 1:num_exps
    exp(e).ensemble = double(ncread(nc(e).obs_epoc, 'observations'));
    exp(e).obs_ind  = -1 * double(ncread(nc(e).obs_epoc, 'obs_type'));
    exp(e).O_time   = double(ncread(nc(e).obs_epoc, 'time')) + origin;
    
    % inflation files
    if inf_flav_n(e, 1) > 0 || inf_flav_n(e, 2) > 0
        if inf_flav_n(e, 1) > 0 && inf_flav_n(e, 2) > 0 
            exp(e).pr.infm.x1 = ncread(nc(e).pr_inflate_mean, 'seg_inflow');
            exp(e).pr.infs.x1 = ncread(nc(e).pr_inflate_sd  , 'seg_inflow');
            exp(e).po.infm.x1 = ncread(nc(e).po_inflate_mean, 'seg_inflow');
            exp(e).po.infs.x1 = ncread(nc(e).po_inflate_sd  , 'seg_inflow');
        elseif inf_flav_n(e, 1) > 0 
            exp(e).pr.infm.x1 = ncread(nc(e).pr_inflate_mean, 'seg_inflow');
            exp(e).pr.infs.x1 = ncread(nc(e).pr_inflate_sd  , 'seg_inflow');
        elseif inf_flav_n(e, 2) > 0 
            exp(e).po.infm.x1 = ncread(nc(e).po_inflate_mean, 'seg_inflow');
            exp(e).po.infs.x1 = ncread(nc(e).po_inflate_sd  , 'seg_inflow');
        end

        % Hybrid files
        if hyb_flav_n(e) > 0 
            exp(e).pr.hybm.x1 = ncread(nc(e).hybrid_mean, 'seg_inflow');
            exp(e).pr.hybs.x1 = ncread(nc(e).hybrid_sd  , 'seg_inflow');
        end
    end
end
if diag_ol
    ol.ensemble = double(ncread(ol.obs_epoc, 'observations')); 
    ol.obs_ind  = -1 * double(ncread(ol.obs_epoc, 'obs_type'));
    ol.O_time   = double(ncread(ol.obs_epoc, 'time')) + origin;
end

openloop    = cell( 4, gauges.want.num          );
forecast    = cell(12, gauges.want.num, num_exps);
observation = cell( 9, gauges.want.num, num_exps);
analysis    = cell( 9, gauges.want.num, num_exps);

% prior+posterior ensemble + 2 mean and 2 sd copies
% + 1 obs copy and 1 obs_var copy
Ne = zeros(1, num_exps);
if diag_ol, ol_Ne = ol.ens_size * 2 + 6; end

for e = 1:num_exps
    Ne(e) = exp(e).ens_size * 2 + 6;
end

loc.yo = 1;
loc.mf = 2;
loc.ma = 3;
loc.sf = 4;
loc.sa = 5;

bad_gages = zeros(1, gauges.want.num);

for i = 1:gauges.want.num

    k = gauges.want.IND(i);
    
    for e = 1:num_exps
        loc.ef = 6:2:Ne(e)-1;
        loc.ea = 7:2:Ne(e);
        loc.vo = Ne(e);
    
        exp_obs = k == exp(e).obs_ind;

        fprintf('exp: %2d, obs no: %3d,  dart-index: %6d,  POI gauge ID: %10d\n', e, i, k, gauges.want.OID(i))

        tmp.obs_val    = exp(e).ensemble(loc.yo, exp_obs);
        tmp.ens_mean_f = exp(e).ensemble(loc.mf, exp_obs);
        tmp.ens_sd_f   = exp(e).ensemble(loc.sf, exp_obs);
        tmp.ensemble_f = exp(e).ensemble(loc.ef, exp_obs);
        tmp.obs_var    = exp(e).ensemble(loc.vo, exp_obs); 
        tmp.ens_mean_a = exp(e).ensemble(loc.ma, exp_obs);
        tmp.ens_sd_a   = exp(e).ensemble(loc.sa, exp_obs);
        tmp.ensemble_a = exp(e).ensemble(loc.ea, exp_obs);
        
        if diag_ol
            loc.ol_ef = 6:2:ol_Ne-1;

            % Just in case, open loop has more gauges
            ol_obs  = k == ol.obs_ind;
            ol_time = ol.O_time(ol_obs);

            tmp.ens_mean_ol = ol.ensemble(loc.mf   , ol_obs); 
            tmp.ens_sd_ol   = ol.ensemble(loc.sf   , ol_obs); 
            tmp.ensemble_ol = ol.ensemble(loc.ol_ef, ol_obs);

            if isempty(ol_time)
                bad_gages(i) = 1;
                continue 
            end
        end  

        exp_time = exp(e).O_time(exp_obs); 
        if isempty(exp_time) 
            bad_gages(i) = 1;
            continue 
        end

        % time management:
        [el, ~] = ismember(current, exp_time); 
        
        % check if there is a time-offset issue or 
        % just small number of data from this gauge
        gauge_data_in_time = sum(el);
        if gauge_data_in_time >= 0
            % account for any time shift
            td = current(1) - exp_time(1);
            exp_time = exp_time + td;
            [~, et] = ismember(current, exp_time); et(et == 0) = NaN;
            
            if diag_ol 
                ol_time = ol_time + td; 
                [~, ot] = ismember(current,  ol_time); ot(ot == 0) = NaN;   
            end
        else
            % Something is off with the time between 
            % obs_diag and obs_epoch
            bad_gages(i) = 1;
            continue
        end
    
        obs_val    = zeros(1, Nt);
        ens_mean_f = zeros(1, Nt);
        ens_sd_f   = zeros(1, Nt);
        ensemble_f = NaN(exp(e).ens_size, Nt);
        obs_var    = zeros(1, Nt);
        ens_mean_a = zeros(1, Nt);
        ens_sd_a   = zeros(1, Nt);
        ensemble_a = NaN(exp(e).ens_size, Nt);

        % exp-loop
        for j = 1:Nt
            if j < Nt && ~isnan(et(j)) && ~isnan(et(j+1))
                obs_val(j)       = mean(tmp.obs_val   (:, et(j):et(j+1)), 2);
                ens_mean_f(j)    = mean(tmp.ens_mean_f(:, et(j):et(j+1)), 2);
                ens_sd_f(j)      = mean(tmp.ens_sd_f  (:, et(j):et(j+1)), 2);
                ensemble_f(:, j) = mean(tmp.ensemble_f(:, et(j):et(j+1)), 2);
                obs_var(j)       = mean(tmp.obs_var   (:, et(j):et(j+1)), 2);
                ens_mean_a(j)    = mean(tmp.ens_mean_a(:, et(j):et(j+1)), 2);
                ens_sd_a(j)      = mean(tmp.ens_sd_a  (:, et(j):et(j+1)), 2);
                ensemble_a(:, j) = mean(tmp.ensemble_a(:, et(j):et(j+1)), 2);

            elseif ~isnan(et(j))
                obs_val(j)       = tmp.obs_val   (:, et(j));
                ens_mean_f(j)    = tmp.ens_mean_f(:, et(j));
                ens_sd_f(j)      = tmp.ens_sd_f  (:, et(j));
                ensemble_f(:, j) = tmp.ensemble_f(:, et(j));
                obs_var(j)       = tmp.obs_var   (:, et(j));
                ens_mean_a(j)    = tmp.ens_mean_a(:, et(j));
                ens_sd_a(j)      = tmp.ens_sd_a  (:, et(j));
                ensemble_a(:, j) = tmp.ensemble_a(:, et(j));

            elseif j < Nt && ~isnan(et(j+1))
                obs_val(j)       = tmp.obs_val   (:, et(j+1));
                ens_mean_f(j)    = tmp.ens_mean_f(:, et(j+1));
                ens_sd_f(j)      = tmp.ens_sd_f  (:, et(j+1));
                ensemble_f(:, j) = tmp.ensemble_f(:, et(j+1));
                obs_var(j)       = tmp.obs_var   (:, et(j+1));
                ens_mean_a(j)    = tmp.ens_mean_a(:, et(j+1));
                ens_sd_a(j)      = tmp.ens_sd_a  (:, et(j+1));
                ensemble_a(:, j) = tmp.ensemble_a(:, et(j+1));

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

        % open loop
        if diag_ol
            ens_mean_ol = zeros(1, Nt); 
            ens_sd_ol   = zeros(1, Nt); 
            ensemble_ol = NaN(ol.ens_size, Nt);
            
            for j = 1:Nt
                if j < Nt && ~isnan(ot(j)) && ~isnan(ot(j+1))
                    ens_mean_ol(j)    = mean(tmp.ens_mean_ol(:, ot(j):ot(j+1)), 2); 
                    ens_sd_ol(j)      = mean(tmp.ens_sd_ol  (:, ot(j):ot(j+1)), 2); 
                    ensemble_ol(:, j) = mean(tmp.ensemble_ol(:, ot(j):ot(j+1)), 2);
    
                elseif ~isnan(ot(j))
                    ens_mean_ol(j)    = tmp.ens_mean_ol(:, ot(j));
                    ens_sd_ol(j)      = tmp.ens_sd_ol  (:, ot(j));
                    ensemble_ol(:, j) = tmp.ensemble_ol(:, ot(j));
                        
    
                elseif j < Nt && ~isnan(ot(j+1))
                    ens_mean_ol(j)    = tmp.ens_mean_ol(:, ot(j+1));
                    ens_sd_ol(j)      = tmp.ens_sd_ol  (:, ot(j+1));
                    ensemble_ol(:, j) = tmp.ensemble_ol(:, ot(j+1));
    
                else
                    ens_mean_ol(j)    = NaN; 
                    ens_sd_ol(j)      = NaN; 
                    ensemble_ol(:, j) = NaN;
                end
            end
        end
        clear tmp
        
        varname_f = strcat(state_tag, sprintf('%06d', gauges.want.IND(i)), '_guess');
        varname_a = strcat(state_tag, sprintf('%06d', gauges.want.IND(i)), '_analy');
    
        % Manage open loop copies
        if diag_ol, rmse_ol = abs(ens_mean_ol - obs_val); end
    
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
        if diag_ol
            openloop(:, i)   = { ens_mean_ol    ; ...
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
gage_list = 1:gauges.want.num;
gage_list = gage_list(~bad_gages);
num_gauge = length(gage_list);

if disp_res
    pos1 = 200;
    pos2 = 600;
    
    figL = 750;
    figW = 450;
    
    if num_exps == 2
        figL = figL * 2;
    elseif num_exps > 2
        figL = figL * 2;
        figW = figW * 2;
    end 
    rows = ceil(num_exps/2);
    cols = 1 + (num_exps > 1);

    % All gauges in a single region are appended together.
    % If the file exists, then remove it
    gages_fig_name = strcat(figs_dir, 'all_gauges.pdf');
    if isfile(gages_fig_name), delete(gages_fig_name); end

    disp(' ')
    disp('PLOTTING:')
    disp('=========')

    figure('pos', [pos1, pos2, figL, figW])

    % TIME SERIES EVOLUTION:
    ii = 0;
    for o = gage_list
        ii = ii + 1;

        fprintf('Gauge[%3d]: %10d; %3d/%3d\n', o, gauges.want.OID(o), ii, num_gauge)
        
        for e = 1:num_exps
            subplot(rows, cols, e)

            en_f = plot(current, forecast{8, o, e}   , '-', 'Color', gY); hold on
            en_a = plot(current, analysis{7, o, e}   , '-', 'Color', lB); grid on 
            
            ob_a = plot(current, observation{7, o, e}, '*', 'Color', gR); 
            ob_r = plot(current, observation{8, o, e}, '*', 'Color', rD); 
            
            if diag_ol, op   = plot(current, openloop{1, o}      , '-', 'Color', oR, 'LineWidth', 3); end
            mf   = plot(current, forecast{6, o, e}   , '-', 'Color', bK, 'LineWidth', 3); 
            ma   = plot(current, analysis{6, o, e}   , '-', 'Color', bL, 'LineWidth', 3); 

            if diag_ol, so   = plot(current, openloop{3, o}      , '--', 'Color', oR, 'LineWidth', 1); end
            sf   = plot(current, forecast{4, o, e}   , '--', 'Color', bK, 'LineWidth', 1); 
            sa   = plot(current, analysis{4, o, e}   , '--', 'Color', bL, 'LineWidth', 1); 
            
            limsy = get(gca, 'YLim');
            
            strY = 'Sreamflow (cfs)'; 

            set(gca, 'FontSize', 16, 'XLim', [xticks(1), xticks(end)], 'XTick', xticks, 'XTickLabel', xtickslabel, ...
                     'Ylim', [0 limsy(2)])
            ylabel(strY, 'FontSize', 18)

            if mean(observation{9, o, e}, 'omitnan') < obserr_L
                obs_status = 'Used Obs (Assimilated)';
            else
                obs_status = 'Used Obs (Evaluated only)';
            end

            if diag_ol
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
            if num_exps > 1, title(L, exp_name(e), 'FontSize', 14); end

            str1 = [ 'Segment Index: ', num2str(gauges.want.IND(o)), ' | Gauge ID: ', num2str(gauges.want.OID(o)) ];
            title(str1, 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none')
        end
        
        % save the whole regions
        exportgraphics(gcf, gages_fig_name, 'Append', true)
    
        % separate figures for each gauge
        gag_name = strcat(figs_dir, num2str(gauges.want.OID(o)), form);
        exportgraphics(gcf, gag_name)
        clf
    end
    close
end 

% clean-up!
openloop    = openloop(:, ~bad_gages      );
forecast    = forecast(:, ~bad_gages,    :);
observation = observation(:, ~bad_gages, :);
analysis    = analysis(:, ~bad_gages,    :);

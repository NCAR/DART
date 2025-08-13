clear 
close all
clc

Wong_cols = [ ...
      0,   0,   0; ...
    230, 159,   0; ...
     86, 180, 233; ...
      0, 158, 115; ...
    240, 228,  66; ...
      0, 114, 178; ...
    213,  94,   0; ...
    204, 121, 167; ...
    ] / 255;

source_dir = '/PATH/TO/PWS-DART/EXPERIMENTS/DIRECTORY/'; 
exps_dir   = strcat(source_dir, 'pywatershed/');
dir_ol     = {strcat(exps_dir, 'OL')};
dir_exps   = {strcat(exps_dir, 'DA')};

% Where to save the figures
figs_dir = strcat(exps_dir, 'hydrographs/');

% image format (only for separate gauges and NOT regions)
form = '.png'; %png, jpg, tif, emf, eps

% Run pwsDART diags
[obs_all, openloop_all, forecast_all, analysis_all] = pwsDARTdiags(dir_exps, dir_ol, 0, figs_dir);

% Separate exp name from path
exp_name = string(missing);
for e = 1:length(dir_exps)
    sp_names    = strsplit(dir_exps{e}, '/');
    exp_name(e) = sp_names{end};
end

exp_names_labels = ['OL', exp_name];

%% Compute RMSE-based on percentiles

any_exp = 1;
obs_val = 6;
metrics = 2;
bias_dg = 3;

pr_l = 10; %10th percentile
pr_h = 90; %90th percentile

Ny = size(forecast_all, 2); 
Ne = length(dir_exps);

OBS_pl = zeros(Ny, 1);
OBS_ph = zeros(Ny, 1);

for e = 1:Ne+1
    for k = 1:Ny  
        gage_dat = obs_all{obs_val, k, any_exp};

        OBS_pl(k) = prctile(gage_dat, pr_l);
        OBS_ph(k) = prctile(gage_dat, pr_h);
    
        low_flow = gage_dat <= OBS_pl(k); 
        hig_flow = gage_dat >= OBS_ph(k);
    
        if e == 1
            rmse.f = openloop_all{metrics, k};
            rmse.a = rmse.f;
        else
            rmse.f = forecast_all{metrics, k, e-1};
            rmse.a = analysis_all{metrics, k, e-1};
            bias.f = forecast_all{bias_dg, k, e-1};
            bias.a = analysis_all{bias_dg, k, e-1};
        end
    
        p_RMSE_LOW.f(k, e) = mean(rmse.f(low_flow), 'omitnan');
        p_RMSE_LOW.a(k, e) = mean(rmse.a(low_flow), 'omitnan');

        p_RMSE_HIG.f(k, e) = mean(rmse.f(hig_flow), 'omitnan');
        p_RMSE_HIG.a(k, e) = mean(rmse.a(hig_flow), 'omitnan');

        p_RMSE_ALL.f(k, e) = mean(rmse.f, 'omitnan');
        p_RMSE_ALL.a(k, e) = mean(rmse.a, 'omitnan');

        if e > 1
            p_BIAS_ALL.f(k, e-1) = mean(bias.f, 'omitnan');
            p_BIAS_ALL.a(k, e-1) = mean(bias.a, 'omitnan');
        end
    end
end

%% KGE + NSE data
ole_avg = 1;
ens_avg = 6;

KGE    = zeros(Ny, 2, Ne+1); 
NSE    = zeros(Ny, 2, Ne+1);  
CV_obs = zeros(Ny, Ne+1);
CV_sim = zeros(Ny, Ne+1);
GAMMA  = zeros(Ny, Ne+1);
for e = 1:Ne+1
    for k = 1:Ny  
        Q_obs   = obs_all{obs_val, k, any_exp};
        mu_obs  = mean(Q_obs, 'omitnan');
        sig_obs = std(Q_obs, 'omitnan'); 

        if e == 1
            Q_sim = openloop_all{ole_avg, k};
            Q_upd = Q_sim;
        else
            Q_sim = forecast_all{ens_avg, k, e-1};
            Q_upd = analysis_all{ens_avg, k, e-1};
        end
        mu_sim  = mean(Q_sim, 'omitnan');
        sig_sim = std(Q_sim, 'omitnan');

        mu_s_a  = mean(Q_upd, 'omitnan');
        sig_s_a = std(Q_upd, 'omitnan');

        alpha = sig_sim / sig_obs;
        beta  = mu_sim / mu_obs;

        alpha_a = sig_s_a / sig_obs;
        beta_a  = mu_s_a / mu_obs;

        R = corrcoef(Q_obs, Q_sim, 'Rows', 'complete');
        r = R(1, 2);

        R = corrcoef(Q_obs, Q_upd, 'Rows', 'complete');
        r_a = R(1, 2);

        NSE(k, 1, e) = 1 - sum((Q_sim - Q_obs).^2, 'omitnan') / sum((Q_obs - mu_obs).^2, 'omitnan');
        NSE(k, 2, e) = 1 - sum((Q_upd - Q_obs).^2, 'omitnan') / sum((Q_obs - mu_obs).^2, 'omitnan');

        KGE(k, 1, e) = 1 - sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2);
        KGE(k, 2, e) = 1 - sqrt((r_a-1)^2 + (alpha_a-1)^2 + (beta_a-1)^2);

        CV_obs(k, e) = sig_obs / mu_obs;
        CV_sim(k, e) = sig_sim / mu_sim;
        GAMMA (k, e) = CV_sim(k, e) / CV_obs(k, e);
    end
end


%% RMSE Figs - low+high
figure('pos', [285, 420, 750, 750]) 

offset = 0.15;
bp_wid = 0.3;
bp_sym = '+';  
bp_wis = 1;

subplot('position', [0.1, 0.51, 0.8, 0.4])

vars = p_RMSE_LOW.f; Ns = size(vars, 2);
xmin = min(vars(:));
xmax = max(vars(:));

COLS = brewermap(Ns+1, 'Set1');
COLS(6, :) = Wong_cols(1, :);

boxplot(vars, 'Whisker', bp_wis, 'Widths', bp_wid, 'OutlierSize', 8, ...
    'Colors', COLS, 'Symbol', bp_sym)

set(findobj(gca,'type','line'), 'linew', 2)

set(gca, 'FontSize', 16, 'XLim', [.5, Ns+0.5], 'XTick', 1:Ns, 'XTickLabel', {}, ...
     'YLim', [.01*xmin, 100*xmax], 'YScale', 'log', ...
     'XGrid', 'on', 'YTickLabelRotation', 90)

ylabel('Low Flow Periods', 'FontSize', 18)

title({'HydroDART Summarized Diagnostics'; ...
      ['Prior RMSE, averaged in time and over ' num2str(Ny) ' gauges']}, 'FontSize', 20)

mRMSE_L = mean(vars, 1);

yy = 0.1 * xmin;
if xmin <= 0 
    yy = 2 * xmax;
end
for k = 1:Ns
    text(k-offset, yy, sprintf('%.2f', mRMSE_L(k)), 'FontSize', 14, 'Color', COLS(k, :), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

% high ones
subplot('position', [0.1, 0.1, 0.8, 0.39])

vars = p_RMSE_HIG.f; Ns = size(vars, 2);
xmin = min(vars(:));
xmax = max(vars(:));

COLS = brewermap(Ns+1, 'Set1');
COLS(6, :) = Wong_cols(1, :);

boxplot(vars, 'Whisker', bp_wis, 'Widths', bp_wid, 'OutlierSize', 8, ...
    'Colors', COLS, 'Symbol', bp_sym)

set(findobj(gca,'type','line'), 'linew', 2)

set(gca, 'FontSize', 16, 'XLim', [.5, Ns+0.5], 'XTick', 1:Ns, 'XTickLabel', ...
    exp_names_labels, 'YLim', [.01*xmin, 100*xmax], 'YScale', 'log', ...
     'XGrid', 'on', 'YTickLabelRotation', 90)

ylabel('High Flow Periods', 'FontSize', 18)

mRMSE_H = mean(vars, 1);

yy = 0.1 * xmin;
if xmin <= 0 
    yy = 2 * xmax;
end
for k = 1:Ns
    text(k-offset, yy, sprintf('%.2f', mRMSE_H(k)), 'FontSize', 14, 'Color', COLS(k, :), ...
        'FontSize', 16, 'FontWeight', 'bold');
end
fig_name = strcat(figs_dir, 'RMSE_low_high', form);
exportgraphics(gcf, fig_name)

%% RMSE Figs - all gauges
figure('pos', [285, 420, 600, 750]) 

subplot('position', [0.1, 0.51, 0.8, 0.4])

vars = p_RMSE_ALL.f; Ns = size(vars, 2);
xmin = min(vars(:));
xmax = max(vars(:));

COLS = brewermap(Ne+1, 'Set1');
COLS(6, :) = Wong_cols(1, :);

boxplot(vars, 'Whisker', bp_wis, 'Widths', bp_wid, 'OutlierSize', 8, ...
    'Colors', COLS, 'Symbol', bp_sym)

set(findobj(gca,'type','line'), 'linew', 2)

set(gca, 'FontSize', 16, 'XLim', [.5, Ns+0.5], 'XTick', 1:Ns, 'XTickLabel', {}, ...
     'YLim', [.01*xmin, 100*xmax], 'YScale', 'log', ...
     'XGrid', 'on', 'YTickLabelRotation', 90)

ylabel('Prior Diags', 'FontSize', 18)

title('RMSE for All gauges, averaged over time', 'FontSize', 20)

mRMSE = mean(vars, 1);

yy = 0.3 * xmin;
if xmin <= 0 
    yy = 2 * xmax;
end
for k = 1:Ns
    text(k-offset, yy, sprintf('%.2f', mRMSE(k)), 'FontSize', 16, 'Color', COLS(k, :), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

subplot('position', [0.1, 0.1, 0.8, 0.39])

vars = p_RMSE_ALL.a; Ns = size(vars, 2);
xmin = min(vars(:));
xmax = max(vars(:));

COLS = brewermap(Ne+1, 'Set1');
COLS(6, :) = Wong_cols(1, :);

boxplot(vars, 'Whisker', bp_wis, 'Widths', bp_wid, 'OutlierSize', 8, ...
    'Colors', COLS, 'Symbol', bp_sym)

set(findobj(gca,'type','line'), 'linew', 2)

set(gca, 'FontSize', 16, 'XLim', [.5, Ns+0.5], 'XTick', 1:Ns, 'XTickLabel', ...
    exp_names_labels, 'YLim', [.01*xmin, 100*xmax], 'YScale', 'log', ...
     'XGrid', 'on', 'YTickLabelRotation', 90)

ylabel('Posterior Diags', 'FontSize', 18)

mRMSE = mean(vars, 1);

yy = 0.3 * xmin;
if xmin <= 0 
    yy = 2 * xmax;
end
for k = 1:Ns
    text(k-offset, yy, sprintf('%.2f', mRMSE(k)), 'FontSize', 16, 'Color', COLS(k, :), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

fig_name = strcat(figs_dir, 'RMSE_all', form);
exportgraphics(gcf, fig_name)

%% KGE boxplots
figure('pos', [285, 420, 600, 750]) 

subplot('position', [0.1, 0.51, 0.8, 0.4])

offset = 0.20;

vars = squeeze(KGE(:, 2, :)); Ns = size(vars, 2);

COLS = brewermap(Ne+1, 'Set1');
COLS(6, :) = Wong_cols(1, :);

boxplot(vars, 'Whisker', bp_wis, 'Widths', bp_wid, 'OutlierSize', 8, ...
    'Colors', COLS, 'Symbol', bp_sym)

set(findobj(gca,'type','line'), 'linew', 2)
yl = max(-40, 1.5*min(vars(:)));

set(gca, 'FontSize', 16, 'XLim', [.5, Ns+0.5], 'XTick', 1:Ns, 'XTickLabel', {}, ...
    'YLim', [yl, 5], 'YTick', [-30, -20, -10, -5, 0], 'YTickLabelRotation', 90)

title('Summary hydrologic efficiency for all gauges', 'FontSize', 20)
ylabel('KGE', 'FontSize', 18)

mKGE = median(vars);

yy = 3;
for k = 1:Ns
    text(k-offset, yy, sprintf('%.2f', mKGE(k)), 'FontSize', 16, 'Color', COLS(k, :), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

subplot('position', [0.1, 0.1, 0.8, 0.39])

vars = squeeze(NSE(:, 2, :)); Ns = size(vars, 2);

xmin = min(vars(:));
xmax = max(vars(:));

COLS = brewermap(Ne+1, 'Set1');
COLS(6, :) = Wong_cols(1, :);

boxplot(vars, 'Whisker', bp_wis, 'Widths', bp_wid, 'OutlierSize', 8, ...
    'Colors', COLS, 'Symbol', bp_sym)

set(findobj(gca,'type','line'), 'linew', 2)
yl = max(-40, 1.5*min(vars(:)));

set(gca, 'FontSize', 16, 'XLim', [.5, Ns+0.5], 'XTick', 1:Ns, 'XTickLabel', ...
    exp_names_labels, 'YLim', [yl, 5], 'YTick', [-30, -20, -10, -5, 0], 'YTickLabelRotation', 90)

ylabel('NSE', 'FontSize', 18)

mKGE = median(vars);

yy = 3;
for k = 1:Ns
    text(k-offset, yy, sprintf('%.2f', mKGE(k)), 'FontSize', 16, 'Color', COLS(k, :), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

fig_name = strcat(figs_dir, 'KGE_NSE_all', form);
exportgraphics(gcf, fig_name)



%% Collect data in file

% 4 files: low, high, prior, post (RMSE changes per file)
% +4 more for KGE and NSE of priors and posteriors
% gauge id RMSE-ol RMSE-enkf RMSE-h1 RMSE-h2 ...

obs_id = 2;

% gage IDs:
gages = zeros(Ny, 1);
for k = 1:Ny  
    gages(k) = obs_all{obs_id, k, any_exp};
end

% files
f.prior = strcat(figs_dir, 'prior_diags.txt');
f.poste = strcat(figs_dir, 'poste_diags.txt');
f.highf = strcat(figs_dir, 'highf_diags.txt');
f.lowfl = strcat(figs_dir, 'lowfl_diags.txt');

f.kge_f = strcat(figs_dir, 'KGE_f_diags.txt');
f.kge_a = strcat(figs_dir, 'KGE_a_diags.txt');
f.nse_f = strcat(figs_dir, 'NSE_f_diags.txt');
f.nse_a = strcat(figs_dir, 'NSE_a_diags.txt');

% data
d.prior = [gages, p_RMSE_ALL.f]';
d.poste = [gages, p_RMSE_ALL.a]';
d.highf = [gages, p_RMSE_HIG.f]';
d.lowfl = [gages, p_RMSE_LOW.f]';

d.kge_f = [gages, squeeze(KGE(:, 1, :))]';
d.kge_a = [gages, squeeze(KGE(:, 2, :))]';
d.nse_f = [gages, squeeze(NSE(:, 1, :))]';
d.nse_a = [gages, squeeze(NSE(:, 2, :))]';

% write to file
s = repmat(' %12.5f ', 1, Ne+1);
s = ['%10d', s, '\n'];

o = fopen(f.prior, 'w');
fprintf(o, s, d.prior);
fclose(o);

o = fopen(f.poste, 'w');
fprintf(o, s, d.poste);
fclose(o);

o = fopen(f.highf, 'w');
fprintf(o, s, d.highf);
fclose(o);

o = fopen(f.lowfl, 'w');
fprintf(o, s, d.lowfl);
fclose(o);

o = fopen(f.kge_f, 'w');
fprintf(o, s, d.kge_f);
fclose(o);

o = fopen(f.kge_a, 'w');
fprintf(o, s, d.kge_a);
fclose(o);

o = fopen(f.nse_f, 'w');
fprintf(o, s, d.nse_f);
fclose(o);

o = fopen(f.nse_a, 'w');
fprintf(o, s, d.nse_a);
fclose(o);

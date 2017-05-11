function plot_reg_factor(fname)
%% plot_reg_factor

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (nargin > 1)
    error('Must supply at most 1 filename.')
elseif (nargin ==1)
    % nothing do do
else
    fname = input('Input name of regression diagnostics file:\n<cr> for reg_diagnostics\n','s');
    if isempty(fname)
        fname = 'reg_diagnostics';
    end
end

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

obs_index = input('Input observation index: \n');

disp('loading usually takes a while ... please be patient ...')
load(fname)
disp('done loading.')

% write(diag_unit, 22) days, secs, obs_index, state_index, comp_reg_factor

t     = reg_diagnostics(:,1) + reg_diagnostics(:,2)/86400.0;
dt    = diff(t);
tinds = 1 + find(dt > max(dt)/2);    % indices of the new epoch.
num_times = length(tinds) + 1;       % Number of time steps.

  obs_indices = reg_diagnostics(1:tinds(1)-1,3);
state_indices = reg_diagnostics(1:tinds(1)-1,4);

dobs        = diff(obs_indices);
dinds       = 1 + find(dobs > max(dobs)/2);
num_obs     = length(dinds) + 1;

num_state = size(reg_diagnostics,1)/(num_obs*num_times);

  obs_ind_mat = reshape(  obs_indices,num_state,num_obs);
state_ind_mat = reshape(state_indices,num_state,num_obs);

unique_obs    =   obs_ind_mat(1,:);   % matrix of obs   indices
unique_states = state_ind_mat(:,1);   % matrix of state indices

% This banks on the fact the state index varies fastest, then obs, then time.

plotdat.data = reshape(reg_diagnostics(:,5),[num_state num_obs num_times]);

subset = squeeze(plotdat.data(:,obs_index,:));

plotdat.shapes = {'states','observations','times'};
plotdat.x      = unique_states;
plotdat.mean   = mean(subset,2);
plotdat.median = median(subset,2);
plot(plotdat.x,plotdat.mean,'k-',plotdat.x, plotdat.median,'r-')
title('regression factor diagnostics')
xlabel('state variable index')

legend('mean','median')
legend boxoff


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

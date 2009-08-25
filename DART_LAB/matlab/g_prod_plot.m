function [prior_mean, prior_sd, obs_mean, obs_err_sd, is_err] = g_prod_plot(h)
% Updates the plot of the prior and observation gaussians

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2009, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Successful return as defaul  t
is_err = false;

% Default failed returns for other quantities
prior_mean = 0; prior_sd   = -1;
obs_mean   = 0; obs_err_sd = -1;

h_prior_mean = get(h.edit1);
prior_mean = str2double(h_prior_mean.String);
% The mean must be a number
if(isnan(prior_mean)) 
   error_banner(h, 'Prior Mean must be a number');
   is_err = true;
   return
end

h_prior_sd = get(h.edit2);
prior_sd = str2double(h_prior_sd.String);

% Prior sd must be a number
if(isnan(prior_sd))
   error_banner(h, 'Prior SD must be a postive number');
   is_err = true;
   return
end

% Prior sd must also be positive 
if(prior_sd <= 0)
   error_banner(h, 'Prior SD must be positive')
   is_err = true;
   return
end
   
hold off
prior_handle = plot_gaussian(prior_mean, prior_sd, 1);
set(prior_handle, 'Color', [0 0.73 0], 'LineWidth', 2);
hold on

h_obs_mean = get(h.edit3);
obs_mean = str2double(h_obs_mean.String);

% Obs value must be a number
if(isnan(obs_mean)) 
   error_banner(h, 'Obs value must be a number');
   is_err = true;
   return
end

h_obs_err_sd = get(h.edit4);
obs_err_sd = str2double(h_obs_err_sd.String);

% Obs error sd must be a positive number
if(isnan(obs_err_sd))
   error_banner(h, 'Obs Error SD must be a positive number');
   is_err = true;
   return
end

if(obs_err_sd <= 0)
   error_banner(h, 'Obs Error SD must be positive');
   is_err = true;
   return
end

obs_handle = plot_gaussian(obs_mean, obs_err_sd, 1);
set(obs_handle, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

% Put on a legend
legend('Prior', 'Obs. Likelihood');

end

%---------------------------------------------------------

% Internal function to write error banner
function error_banner(h, message)

   hold off
   x= [1 2];
   plot(x, 'Visible', 'off')
   h_fig = get(h.figure1);
   x_limits = get(h_fig.CurrentAxes, 'Xlim');
   y_limits = get(h_fig.CurrentAxes, 'Ylim');
   text(x_limits(1) * 7/8 + x_limits(2) / 8, mean(y_limits), ...
      message, 'fontsize', 16, 'Color', 'r');

   % While we're at it, clear out the values for posterior, too
   set(h.text7, 'String', 'Posterior Mean = ');
   set(h.text8, 'String', 'Posterior SD = ');
   set(h.text9, 'String', 'Weight = ');
   return; 

end

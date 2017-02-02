function [prior_mean, prior_sd, obs_mean, obs_err_sd, is_err] = g_prod_plot(h)
%% g_prod_plot Updates the plot of the prior and observation gaussians

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Successful return as default
is_err = false;

atts = stylesheet;   % get the uniform colors, etc.

FontSize = get(h.graph,'FontSize');

% Default failed returns for other quantities
prior_mean = 0; prior_sd   = -1;
obs_mean   = 0; obs_err_sd = -1;

h_prior_mean = get(h.ui_edit_prior_mean);
prior_mean = str2double(h_prior_mean.String);
% The mean must be a number
if(isnan(prior_mean))
   error_banner(h, 'Prior Mean must be a number');
   is_err = true;
   return
end

h_prior_sd = get(h.ui_edit_prior_sd);
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
set(prior_handle, 'Color', atts.green, 'LineWidth', 2);
hold on

h_obs_mean = get(h.ui_edit_observation);
obs_mean = str2double(h_obs_mean.String);

% Obs value must be a number
if(isnan(obs_mean))
   error_banner(h, 'Obs value must be a number');
   is_err = true;
   return
end

h_obs_err_sd = get(h.ui_edit_obs_error_sd);
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
set(obs_handle, 'Color', atts.red, 'LineStyle', '--', 'LineWidth', 2);

% Put on a title and a legend
title('gaussian_product','Interpreter','none')
legend('Prior', 'Obs. Likelihood', 'Location', 'NorthWest');
set(gca, 'FontSize', atts.fontsize);

end

%---------------------------------------------------------

% Internal function to write error banner
function error_banner(h, message)

   atts = stylesheet;   % get the uniform colors, etc.

   hold off
   x= [1 2];
   plot(x, 'Visible', 'off')
   h_fig = get(h.figure1);
   x_limits = get(h_fig.CurrentAxes, 'Xlim');
   y_limits = get(h_fig.CurrentAxes, 'Ylim');
   text(x_limits(1) * 7/8 + x_limits(2) / 8, mean(y_limits), ...
      message, 'FontSize', atts.fontsize + 2, 'Color', 'r');

   % the values for the posterior get reset by gaussian_product.m 
   return;

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

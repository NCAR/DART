%% DART:script18 Looking at inflation to deal with errors

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Begin by setting random seed to a nice controlled
% initial value that gives nice plots
randn('state', 0);

x_prior = [-2.2 -1.3 -0.4 1.2 2.2];
% Adjust this to have an exact mean and variance consistent with desired distribution
x_prior = (x_prior - mean(x_prior)) * 0.6 / std(x_prior) - 2.0

% Plot the prior distribution
y_prior = [0.02, 0.02, 0.02, 0.02, 0.02];
h_plot = plot(x_prior, y_prior, 'g*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
set(h_plot, 'color', [0 0.73 0]);

axis([-4 4 0 0.8]);
set(gca, 'fontsize', 24);
p_text = text(-2.0, 0.1, 'Prior Ensemble', 'fontsize', 24);
% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 10.5 5.0]);
set(gca, 'linewidth', 2);
grid on;
% Set the shape of the plot box
pbaspect([2.4 1 1]);
ylabel('Probability');
hold on;

% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');

% Plot a gaussian fit to ensemble
set(p_text, 'visible', 'off');
x = -5:0.01:5;
prior = normpdf(x, -2.0, 0.6);
obs = normpdf(x, 2.0, 0.6);
product = prior .* obs;
h_prior = plot(x, prior, 'g', 'linewidth', 3)
set(h_prior, 'color', [0 0.73 0]);
text(-3.6, 0.7, 'Prior PDF', 'fontsize', 24);

% Overlay the S.D. of distribution
sdx = [0, 0.6] - 2.0;
sdy = [0.2, 0.2];
h_sd = plot(sdx, sdy, 'linewidth', 3);
set(h_sd, 'color', [0 0.73 0]);
h_sd_label = text(-2.6, 0.2, 'S.D.', 'fontsize', 24);


plot(x, obs, 'r', 'linewidth', 3)
text(0.2, 0.7, 'Obs. Likelihood', 'fontsize', 24);

% Overlay the S.D. of distribution
sdx = [-0.6, 0] + 2.0;
sdy = [0.2, 0.2];
hhh = plot(sdx, sdy, 'r', 'linewidth', 3);
text(2.1, 0.2, 'S.D.', 'fontsize', 24);

pause;
print -depsc s18f01.eps;

% Now add on the EXPECTED difference given same mean
expected = sqrt(0.6^2 + 0.6^2);
exx = [0, expected] -2.0;
exy = [0.3, 0.3];
hex = plot(exx, exy, 'b', 'linewidth', 3);
h_exp_label = text(-1.1, 0.3, 'Expected Separation', 'fontsize', 24);

pause
print -depsc s18f02.eps;


% Now put on the actual separation with Expected SDs overlapped
ax = [0, 4] - 2.0;
ay = [0.4, 0.4];
ha = plot(ax, ay, 'k', 'linewidth', 3);
% Mark ticks for SDs
for i = 1:6
   tx = [0, 0] - 2.0 + (i - 1) * expected;
   ty = [0.38, 0.42];
   ht(i) = plot(tx, ty, 'b', 'linewidth', 3);
end
diff_sds = 4.0 / expected
sep_label = ['Actual ', num2str(4 / expected), ' SDs'];
h_sep_label = text(-1.2, 0.47, sep_label, 'fontsize', 24);

pause
print -depsc s18f03.eps;

% Now, discuss variance inflation
% First, turn off the old expected and actual for now
%set(h_sd, 'visible', 'off');
%set(h_sd_label, 'visible', 'off');
set(hex, 'visible', 'off');
set(ha, 'visible', 'off');
set(h_sep_label, 'visible', 'off');
set(h_exp_label, 'visible', 'off');
for i = 1:6
   set(ht(i), 'visible', 'off');
end

% Zoom in on the prior portion of the plot

axis([-4 0 0 0.8]);
pause
print -depsc s18f04.eps;


x_prior_inf = 1.5 * (x_prior + 2.0) -2.0;
% Plot the inflated prior distribution
y_prior_inf = y_prior + 0.10;
h_plot = plot(x_prior_inf, y_prior_inf, 'g*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
set(h_plot, 'color', [0 0.73 0]);

% Draw vectors between these guys
for i = 1:5
   x_diff = [x_prior(i), x_prior_inf(i)]; 
   y_diff = [y_prior(i), y_prior_inf(i)]; 
   h_diff(i) = plot(x_diff, y_diff, 'linewidth', 2);
   set(h_diff(i), 'color', [0 0.73 0]);
end

h_inf_text = text(-1.5, 0.7, 'Inflate SD by 1.5', 'fontsize', 24);
h_inf_text2 = text(-1.5, 0.60, 'Variance by 1.5^2', 'fontsize', 24);

pause
print -depsc s18f05.eps;

% Change the original distribution to dashed and plot the new distribution solid
set(h_prior, 'linestyle', '--', 'linewidth', 2);

% Plot a gaussian fit to ensemble
set(p_text, 'visible', 'off');
x = -5:0.01:5;
prior_inf = normpdf(x, -2.0, 0.6 * 1.5);
h_prior_inf = plot(x, prior_inf, 'g', 'linewidth', 3)
set(h_prior_inf, 'color', [0 0.73 0]);

% Set the old SD to dashed
set(h_sd, 'linestyle', '--');
set(h_sd_label, 'visible', 'off');

% Overlay the S.D. of the inflated distribution
sdx = [0, 1.5 * 0.6] - 2.0;
sdy = [0.25, 0.25];
h_sd_inf = plot(sdx, sdy, 'linewidth', 3);
set(h_sd_inf, 'color', [0 0.73 0]);
h_sd_label = text(-3.0, 0.25, 'Inflated S.D.', 'fontsize', 24);
pause
print -depsc s18f06.eps;

% Remove inflate labels
set(h_inf_text, 'visible', 'off');
set(h_inf_text2, 'visible', 'off');

% Now go back to original display
axis([-4 4 0 0.8]);
% Turn off the old S.D. label and move the new SD label to it's place?
set(h_sd, 'visible', 'off');
set(h_sd_label, 'position', [-3.9 0.25 0]);

% Now add on the EXPECTED difference given same mean
expected = sqrt((1.5 *0.6)^2 + 0.6^2);
exx = [0, expected] -2.0;
exy = [0.3, 0.3];
hex = plot(exx, exy, 'b', 'linewidth', 3);
h_exp_label = text(-0.8, 0.3, 'Expected Separation', 'fontsize', 24);

% Now put on the actual separation with Expected SDs overlapped
ax = [0, 4] - 2.0;
ay = [0.4, 0.4];
ha = plot(ax, ay, 'k', 'linewidth', 3);
% Mark ticks for SDs
for i = 1:5
   tx = [0, 0] - 2.0 + (i - 1) * expected;
   ty = [0.38, 0.42];
   ht(i) = plot(tx, ty, 'b', 'linewidth', 3);
end
diff_sds = 4.0 / expected
sep_label = ['Actual ', num2str(4 / expected), ' SDs'];
h_sep_label = text(-1.2, 0.47, sep_label, 'fontsize', 24);

pause
print -depsc s18f07.eps;


% Classical perturbed obs ensemble Kalman filter

% Data Assimilation Research Testbed -- DART
% Copyright 2004-2007, Data Assimilation Research Section
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

x_prior = [-2.2 -1.8 -1.4 1.8 2.2];
% Adjust this to have an exact mean and variance consistent with desired distribution
x_prior = (x_prior - mean(x_prior)) * 1.2 / std(x_prior) - 1.0

% Plot the prior distribution
y_prior = [0.05, 0.05, 0.05, 0.05, 0.05];
hold on;
for i = 1:5
   h_prior_plot(i) = plot(x_prior(i), y_prior(i), 'g*');
   set(h_prior_plot(i), 'markersize', 18);
   set(h_prior_plot(i), 'linewidth', 3);
end

axis([-4 4 0 0.6]);
set(gca, 'fontsize', 24);
p_text = text(-3.8, 0.1, 'Prior Ensemble', 'fontsize', 24);
% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 11 5.0]);
set(gca, 'linewidth', 2);
grid on;
% Set the shape of the plot box
pbaspect([2.4 1 1]);
ylabel('Probability');



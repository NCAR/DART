%% Do some demo plots

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Do ensemble scatter as function of time for a given variable
figure(1);
hold on;
ens_mean = sum(posterior_state(:, :, 1), 2) / copies_per_time;
for i = 1:copies_per_time,
   plot(posterior_state(:, i, 1));
end
plot(ens_mean, 'r');
plot(true_state(:, 1, 1), 'g');



% Do a plot of the ensemble of states at a given time
figure(2);
hold on;
ens_mean = sum(posterior_state(10, :, :), 2) / copies_per_time;
plot(location, squeeze(true_state(10, 1, :)), 'g');

% Plot the ensemble mean, too
plot(location, squeeze(ens_mean), 'r');

for i = 1:copies_per_time,
   plot(location, squeeze(posterior_state(10, i, :)));
end
xlabel('Location');
legend('Ensemble mean', 'Truth', 'Ensembles');
title 'Spatial representation of state at fixed time';



% Look at posterior, prior and truth in mean

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

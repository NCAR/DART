%% option7 Tutorial basics for the rank histogram filter (need better name)
% all input is prompted for ... no functional form of this script.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$
%
% the matlab statistics toolbox is required to run this script.

% Begin by generating a sample draw from a gaussian
close all; clear;

% The observed value and observational error standard deviation
obs_val = 0.5;
obs_err_sd = 1.0;
obs_err_var = obs_err_sd^2;

% Generate an ensemble and sort it
mean = 0;
sd  = 1;
%ens = sort(random('Normal', mean, sd, ens_size, 1));

% Special prior design
ens = [-1.5, -.8, -.3, .3, 1.4];
%ens = [-2.5, -.2, -.12 -.07, 0.0, .05, 0.13, 0.21, .25];
%ens = [-2.1, -2.0, -1.8, 0.3, 0.6, 0.8];
ens_size = size(ens, 2);


% Compute a weight for each of the prior ensemble members from
% definition of gaussian
weight = normpdf(ens, obs_val, obs_err_sd);

% Compute the points that bound all the updated mass boxes; start with ensemble
for i = 1:ens_size
   x(2*i - 1) = ens(i);
end
% Compute the mid-point interior boundaries; halfway between ensembles
for i = 2:2:2*ens_size - 2
   x(i) = (x(i - 1) + x(i + 1)) / 2.0;
end

% Compute the ensemble sample var/sd
ens_var = var(ens);
% The wing boxes are gaussian tails to avoid divergence
% The wing kernels have mean at the outermost ensemble members
% Variance is rather arbitrarily set to total ensemble variance / 2
prior_var = ens_var / 2;
% Total mass in kernel should be 1/(ens_size + 1)
% Only doing half of a gaussian so multiply by 2 / (ens_size + 1)
prior_sd = sqrt(prior_var);

% The flanks are updated by standard gaussian products WITH WEIGHTS
var_ratio = obs_err_var / (prior_var + obs_err_var);
new_var = var_ratio * prior_var;
new_sd = sqrt(new_var);

% Get the left updated mean and weight
new_mean_left =  var_ratio * (ens(1) + prior_var * obs_val / obs_err_var);
new_mean_right = var_ratio * (ens(ens_size) + prior_var * obs_val / obs_err_var);
e = exp(1);
prod_weight_left = e .^ (-0.5 .* (ens(1).^2 ./ prior_var + obs_val.^2 ./ obs_err_var - new_mean_left.^2 ./ new_var)) / sqrt(2 * pi);
prod_weight_right = e .^ (-0.5 .* (ens(ens_size).^2 ./ prior_var + obs_val.^2 ./ obs_err_var - new_mean_right.^2 ./ new_var)) / sqrt(2 * pi);

% Compute the mass in the updated outer regions
mass(1) = normcdf(ens(1), new_mean_left, new_sd) * prod_weight_left * 2.0 / (ens_size + 1);
mass(2*ens_size) = (1.0 -normcdf(ens(ens_size), new_mean_right, new_sd)) * prod_weight_right * 2.0 / (ens_size + 1);

% Compute the mass in the inner half boxes that have ensemble point on left
for i = 2:2:2*ens_size - 2
   mass(i) = (1.0 / (2.0 * (ens_size + 1))) * weight(i/2);
end

% Now right inner half boxes
for i = 3:2:2*ens_size - 1
   mass(i) = (1.0 / (2.0 * (ens_size + 1.0))) * weight(fix(i/2) + 1);
end

% Now normalize the mass in the different bins
mass_sum = sum(mass);
mass = mass / mass_sum;

% Find cumulative mass at each box boundary and middle boundary
% Careful, indexing is different from F90 code due to matlab constraints
cumul_mass(1) = 0.0;
for i = 2:2*ens_size + 1
   cumul_mass(i) = cumul_mass(i - 1) + mass(i - 1);
end





% Get resampled position 
% Need 1/ens_size mass between each
umass = 1.0 / (ens_size + 1);

for i = 1:ens_size
   % If it's in the inner or outer range have to use normal
   if(umass < cumul_mass(2))
      % In the first normal box
      left_weight = (1.0 / mass_sum) * prod_weight_left * (2.0 / (ens_size + 1.0))
      %call weighted_norm_inv(left_weight, new_mean_left, new_sd, umass, new_ens(i))

      % Do the inverse of the normal cdf; normalize
      np = umass / left_weight;
      new_ens(i) = norminv(np, new_mean_left, new_sd);

   elseif(umass > cumul_mass(2*ens_size))
      % In the last normal box; Come in from the outside
      right_weight = (1.0 / mass_sum) * prod_weight_right * (2.0 / (ens_size + 1.0));
      %call weighted_norm_inv(right_weight, new_mean_right, new_sd, 1.0 - umass, new_ens(i))

      % Do the inverse of the normal cdf; normalize
      np = (1.0 - umass) / right_weight;
      new_ens(i) = norminv(np, new_mean_right, new_sd);
      % This one is heading out to the right so flip it
      new_ens(i) = new_mean_right + (new_mean_right - new_ens(i));

   else
      % In one of the inner uniform boxes.
      lowest_box = 1;
      for j = lowest_box: 2 * ens_size - 2
         % Find the box that this mass is in
         if(umass >= cumul_mass(j + 1) & umass <= cumul_mass(j + 2))
            new_ens(i) = x(j) + ((umass - cumul_mass(j + 1)) / (cumul_mass(j+2) - cumul_mass(j + 1))) * (x(j + 1) - x(j));
            % Don't need to search lower boxes again
            lowest_box = j;
            break;
         end 
      end
   end
   % Want equally partitioned mass in update with exception that outermost boxes have half
   umass = umass + 1.0 / (ens_size + 1);
end


%-------- Start plotting section here-----------

% Plot the ensemble on an axis
y(1:ens_size) = -0.15;
h = plot(ens, y, 'g*', 'markersize', 16, 'linewidth', 3);

% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 10.5 5.0]);
set(gca, 'linewidth', 2);
grid on;
% Set the shape of the plot box
pbaspect([2.4 1 1]);

hold on;
set(gca, 'fontsize', 18)

g =text(-2.5, -0.15, 'Prior', 'fontsize', 20);
set(g, 'color', 'g');

ylabel('Probability Density', 'fontsize', 20);
% These setttings are for the base case example only
set(gca, 'YTick', [0, 0.2, 0.4, 0.6]);
axis([-3.5, 3.5, -0.27, 0.6]);

% Put on a black base axis
aa = [-3.5, 3.5];
bb = [0 0];
plot(aa, bb, 'k', 'linewidth', 2);


% Setup the printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print -depsc -tiff option7_1.eps;


% Each box gets 1/(n + 1) of the mass
mass_per_box = 1.0 / (ens_size + 1);

% Plot the density boxes for the prior estimate interior
for i = 1:ens_size - 1
   width(i) = ens(i + 1) - ens(i);
   height(i) = mass_per_box / width(i);
   y = [0, 0, height(i), height(i)];
   xx = [ens(i), ens(i+1), ens(i+1), ens(i)];
   h = fill(xx, y, 'g');
   set(h, 'edgecolor', 'g');
   % Output intermediate box addition
   hhh = figure(1);
   outfile = ['option7_', num2str(i + 1), '.eps'];
   print(hhh, '-depsc', '-tiff', outfile);
end

% Get x-values for left side plot
num_pts = 300;
for i = 1:num_pts
   x_left(i) = ens(1) - ((i - 1) / num_pts) * 4 * prior_sd;
end
y_left = normpdf(x_left, ens(1), prior_sd);
y_left = y_left .* 2 ./ (ens_size + 1);
for i = 1:num_pts - 1
   a(1) = x_left(i);       b(1) = 0;
   a(2) = x_left(i + 1);   b(2) = 0;
   a(3) = a(2);            b(3) = y_left(i+1);
   a(4) = a(1);            b(4) = y_left(i);
   h = fill(a, b, 'g');
   set(h, 'edgecolor', 'g'); 
end

% Now do the right side wing
for i = 1:num_pts
   x_right(i) = ens(ens_size) + ((i - 1) / num_pts) * 4 * prior_sd;
end
y_right = normpdf(x_right, ens(ens_size), prior_sd);
y_right = y_right .* 2 ./ (ens_size + 1);
for i = 1:num_pts - 1
   a(1) = x_right(i);       b(1) = 0;
   a(2) = x_right(i + 1);   b(2) = 0;
   a(3) = a(2);             b(3) = y_right(i+1);
   a(4) = a(1);             b(4) = y_right(i);
   h = fill(a, b, 'g');
   set(h, 'edgecolor', 'g'); 
end

print -depsc -tiff option7_6.eps;


% Plot the weights given from the likelihood
qqq = plot(ens, weight, '*r', 'markersize', 12, 'linewidth', 2);
% Plot a thin line for the whole normal, too
num_pts = 1000;
for i = 1:num_pts
   likex(i) = obs_val - 4 * obs_err_sd + (i / num_pts) * 8 * obs_err_sd;
   likey = normpdf(likex, obs_val, obs_err_sd);
end
plot(likex, likey, 'r');

print -depsc -tiff option7_7.eps;


% Now compute the updated distribution
for i = 1:ens_size
   if(i > 1)
      left_height(i)  = mass(2*i - 1) / (width(i-1) / 2);
      left_x = (ens(i - 1) + ens(i)) / 2;
      right_x = ens(i);
      xx = [right_x, right_x, left_x, left_x];
      y = [0, left_height(i), left_height(i), 0];
      plot(xx,  y, 'b', 'linewidth', 3); 
   end
   if(i < ens_size)
      right_height(i) = mass(2*i) / (width(i) / 2);
      left_x = ens(i);
      right_x = (ens(i) + ens(i + 1)) / 2;
      xx = [right_x, right_x, left_x, left_x];
      y = [0, right_height(i), right_height(i), 0];
      plot(xx,  y, 'b', 'linewidth', 3); 
   end

   % Output intermediate box addition
   hhh = figure(1);
   outfile = ['option7_', num2str(i + 7), '.eps'];
   print(hhh, '-depsc', '-tiff', outfile);
end


% Get x-values for left side plot updated distribution
num_pts = 300;
for i = 1:num_pts
   x_left(i) = ens(1) - ((i - 1) / num_pts) * 4 * prior_sd;
end
left_weight = (1.0 / mass_sum) * prod_weight_left * (2.0 / (ens_size + 1.0));
y_left = normpdf(x_left, new_mean_left, new_sd) * left_weight;
for i = 1:num_pts - 1
   plot(x_left, y_left, 'b', 'linewidth', 3)
end
% Drop a line at the inside edge in case it's neededA
aa = [x_left(num_pts), x_left(num_pts)];
bb = [0, y_left(num_pts)];
plot(aa, bb, 'b', 'linewidth', 3);

% Now do the right side wing
for i = 1:num_pts
   x_right(i) = ens(ens_size) + ((i - 1) / num_pts) * 4 * prior_sd;
end
right_weight = (1.0 / mass_sum) * prod_weight_right * (2.0 / (ens_size + 1.0));
y_right = normpdf(x_right, new_mean_right, new_sd) * right_weight;
for i = 1:num_pts - 1
   plot(x_right, y_right, 'b', 'linewidth', 3);
end
% Drop a line at the inside edge in case it's neededA
aa = [x_right(1), x_right(1)];
bb = [0, y_right(1)];
plot(aa, bb, 'b', 'linewidth', 3);

print -depsc -tiff option7_13.eps

% Plot the new_ens members
y(1:ens_size) = -0.05;
h = plot(new_ens, y, 'b*', 'markersize', 16, 'linewidth', 3);

% Label the update
g =text(-2.5, -0.05, 'Posterior', 'fontsize', 20);
set(g, 'color', 'b');

% Plot lines between them
for i = 1:ens_size
   aa(1) = ens(i);
   aa(2) = new_ens(i);
   bb(1) = -0.15;
   bb(2) = -0.05;
   plot(aa, bb);
end

print -depsc -tiff option7_14.eps


% What would standard EA have done?
prior_mean = sum(ens) / ens_size;
prior_var = var(ens);
prior_sd = sqrt(prior_var);
new_var = 1.0 / (1.0 / prior_var + 1.0 / obs_err_var);
new_sd = sqrt(new_var);
new_mean = new_var * (prior_mean / prior_var + obs_val / obs_err_var);
sd_ratio = new_sd / prior_sd;
eakf_ens = (ens - prior_mean) .* sd_ratio + new_mean;
y(1:ens_size) = -0.25;
plot(eakf_ens, y, 'b*', 'markersize', 16, 'linewidth', 3);

% Plot lines between them
for i = 1:ens_size
   aa(1) = ens(i);
   aa(2) = eakf_ens(i);
   bb(1) = -0.15;
   bb(2) = -0.25;
   plot(aa, bb);
end

% Label the update
g =text(-3.2, -0.25, 'EAKF Posterior', 'fontsize', 20);
set(g, 'color', 'b');

print -depsc -tiff option7_15.eps

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

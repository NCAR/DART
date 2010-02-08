%% DART:script13  A quick look at local regression, etc.

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

% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 6.5 6.5]);
% Set the shape of the plot box
%pbaspect([2.4 1 1]);

% Set up for proper printing characteristics
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual', 'ztickmode', 'manual');

% Take a look at using multiple axes here
r1 = subplot(2, 2, 1);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
ylabel('Unobserved State Variable');
r2 = subplot(2, 2, 2);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
r3 = subplot(2, 2, 4);
set(gca, 'linewidth', 2);
set(gca, 'fontsize', 24);
xlabel('Observed Variable');


% Position the marginal and joint plot boxes
set(r1, 'Position', [0.14, 0.23, 0.10, 0.6550]);
set(r2, 'Position', [0.25, 0.23, 0.6550, 0.6550]);
set(r3, 'Position', [0.25 0.1100 0.6550 0.1000]);

% Turn off the unwanted tick labels
set(r2, 'Xticklabel', []);
set(r2, 'yticklabel', []);
set(r1, 'Xticklabel', []);
set(r3, 'yticklabel', []);


for j = 1:3
   subplot(r2);
   axis([-2 2 0 8]);
   % First look at impacts on a second 'unobserved' variable
   % Use the rnum array to hold something decidely not random
   ens_size = 20;
   if j == 1
      ens_size = 80;
   elseif j == 2
      ens_size = 40;
   elseif j == 3
      ens_size = 20;
   end
   interval = 4 / (ens_size - 1);
   rnum(1, :) = -2.: interval: 2.;
   rnum(2, :) = rnum(1, :) .^2 + 2;
   ens_size = size(rnum, 2);
   % Add in some noise
   noise = 1.0;
   for i = 1:size(rnum, 2)
      rnum(2, i) = rnum(2, i) + normrnd(0, noise);
   end
   



   h_joint_prior = plot(rnum(1, :), rnum(2, :), 'g*');
   set(h_joint_prior, 'markersize', 12);
   set(h_joint_prior, 'linewidth', 2);
   set(h_joint_prior, 'color', [0 0.73 0]);
   grid on;
   set(r2, 'Xticklabel', []);
   set(r2, 'yticklabel', []);
   hold on;
   
   % Plot the y = x^2 curve on which the relation actually lies
   x_limits = get(gca, 'xlim');
   x_pts = x_limits(1) : 0.01: x_limits(2);
   y_pts = x_pts .^ 2 + 2;
   plot(x_pts, y_pts, 'r', 'linewidth', 2);
   
   
   
   % Plot the marginals in the wing plots
   subplot(r3);
   yu(1:ens_size) = 0.1;
   h_prior_marg = plot(rnum(1, :), yu, '*g');
   set(h_prior_marg, 'markersize', 10)
   set(h_prior_marg, 'linewidth', 1.6);
   set(h_prior_marg, 'color', [0 0.73 0]);
   % Use the x axis limits from the main plot
   axis([get(r2, 'xlim')  [0 1]]);
   set(r3, 'yticklabel', []);
   xlabel('Observed Variable');
   
   % Plot the marginals in the wing plots
   subplot(r1);
   xu(1:ens_size) = 0.1;
   h_state_marg = plot(xu, rnum(2, :), '*g');
   set(h_state_marg, 'markersize', 10)
   set(h_state_marg, 'linewidth', 1.6);
   set(h_state_marg, 'color', [0 0.73 0]);
   axis([[0 1] get(r2, 'ylim')]);
   set(r1, 'xticklabel', []);
   ylabel('Unobserved State Variable');
   
   % Now look at fitting regression lines at x 0 for different 
   % numbers of neighbors
   
   
   subplot(r2);
      for i = 1:3
         x_point(1) = (i - 1)*1.25 - 1.75;
         x_point(2) = x_point(1) + 1;
         base_index(i) = 1 + fix((x_point(1) + 2) / interval);
         % Plot the perfect target points
         y_point(2) = x_point(2)^2 + 2;
         h_target = plot(x_point(2), y_point(2), 'b*');
         set(h_target, 'markersize', 12)
         set(h_target, 'linewidth', 2);
         % Plot the marginal prior obs variable as an arrow, too
         subplot(r3);
	 h_arrow = arrow([x_point(1), 0.4], [x_point(2), 0.4]);
        
         subplot(r2);
         % Look at impacts of obs delta 1 at different point
         % Plot best fit line through this mess; just uses the standard deviations?
         bot = base_index(i);
         top = base_index(i) + fix(1 / interval); 
         % compute the sample statistics at a local set of points
         s_correl = corrcoef(rnum(1, bot:top), rnum(2, bot:top));
         s_var(1) = var(rnum(1, bot:top));
         s_var(2) = var(rnum(2, bot:top));
         s_covar = s_correl(1, 2) * sqrt(s_var(1)) * sqrt(s_var(2));
         reg_coef = s_covar / s_var(1);
      
      
         y_point(1) = x_point(1) ^ 2 + 2;
         y_point(2) = y_point(1) + reg_coef;
         if j == 1
            plot(x_point, y_point, 'b', 'linewidth', 3);
            % Fake my own legend
            xu = [-1.75, -1.25]; yu = [7.5, 7.5]; 
            plot(xu, yu, 'b', 'linewidth', 3);
            text(-1.1, 7.5, '80 Member Ensemble', 'fontsize', 24);
         elseif j == 2
            plot(x_point, y_point, 'b--', 'linewidth', 3);
            % Fake my own legend
            xu = [-1.75, -1.25]; yu = [7.0, 7.0]; 
            plot(xu, yu, 'b--', 'linewidth', 3);
            text(-1.1, 7.0, '40 Member Ensemble', 'fontsize', 24);
         elseif j == 3
            plot(x_point, y_point, 'b-.', 'linewidth', 3);
            % Fake my own legend
            xu = [-1.72, -1.25]; yu = [6.5, 6.5]; 
            plot(xu, yu, 'b-.', 'linewidth', 3);
            text(-1.1, 6.5, '20 Member Ensemble', 'fontsize', 24);
         end
      end
   pause
   outfile = ['s13f0', num2str(j), '.eps'];
   print(gcf, '-depsc', outfile);

   clear rnum;
   clear yu;
   clear xu;
   set(h_joint_prior, 'visible', 'off');
   set(h_prior_marg, 'visible', 'off');
end

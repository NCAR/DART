function [xp, yp] = plot_rhf_pdf(x, ens_size, mass, left_mean, left_sd, left_amp, ...
   right_mean, right_sd, right_amp, bounded_left)

% Creates the x and y point vectors to plot a rhf PDF
% Bounds of plotting domain are given by xlow and xhigh
% x is the ensemble values
% mass is the amount of probability mass in each of the ens_size+1 bins
% left_mean, left_sd, and left_amp are the left tail PDF, similar for right

% Use the tail distributions to get 5 standard deviations on the plots
if(bounded_left)
   xlow = 0.0;
else
   xlow = left_mean - 5.0 * left_sd;
end

xhigh = right_mean + 5.0 * right_sd;

% Point spacing
point_del = (xhigh - xlow) / 1000;
% Plot the left tail
xptsb = xlow:point_del:x(1);
yptsb = left_amp * normpdf(xptsb, left_mean, left_sd);

% Do the interior bins
for i = 1:ens_size - 1
   % Do the vertical line on the right
   height(i+1) = mass(i+1) / (x(i + 1) - x(i));

   xpts(3*i - 2) = x(i);
   xpts(3*i - 1) = x(i);
   xpts(3*i)     = x(i + 1);

   if(i == 1) 
      ypts(3*i -2) = left_amp * normpdf(x(1), left_mean, left_sd);
   else
      ypts(3*i -2) = height(i);
   end
   ypts(3*i - 1)   = height(i + 1);
   ypts(3*i)       = height(i+1);
end

% Plot the vertical on right of last interior bin
xptsv(1:2) = x(ens_size);
yptsv(1)   = height(ens_size);
yptsv(2)   = right_amp * normpdf(x(ens_size), right_mean, right_sd);

% Plot the right tail
xptsa = x(ens_size):point_del:xhigh;
yptsa = right_amp * normpdf(xptsa, right_mean, right_sd);

% Concatenate the different pieces
xp = [xptsb xpts xptsv xptsa];
yp = [yptsb ypts yptsv yptsa];

end

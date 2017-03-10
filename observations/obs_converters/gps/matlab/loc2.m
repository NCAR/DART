%% loc2

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

load radi_obs.dat;
err = radi_obs(:, 1);
otype = radi_obs(:, 9);

counter = 0;
for i = 1:size(radi_obs, 1)
  counter = counter + 1;
  xx4(counter) = radi_obs(i, 2)*180.0/3.14159;
  yy4(counter) = radi_obs(i, 3)*180.0/3.14159;
end

 subplot('position', [0.1,0.3,0.8,0.5])
plot(xx4, yy4, 'b.');
axis([90 180 -20 15])            % display only the region seleceted.
ylabel('Latitude (N)', 'fontsize', 12)
xlabel('Longitude (E)', 'fontsize', 12)
title('Radiosonde locations, Jan 6, 2007', 'fontsize', 14)

 subplot('position', [0.0,0.0,0.06,0.06])
 text(0.15, 0.15, 'Fig. 1b')
print -dpsc radi_obs.ps
print -dpng radi_obs.png

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

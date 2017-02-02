%% loc3

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

load gpsloc.20030108;

load comp_obs.20030108;

counter = 0;
for i = 1:size(gpsloc, 1)
  counter = counter + 1;
  xx4(counter) = gpsloc(i, 2)*180.0/3.14159;
  yy4(counter) = gpsloc(i, 3)*180.0/3.14159;
end

counter = 0;
for i = 1:size(comp_obs, 1)
  counter = counter + 1;
  xx4_comp_obs(counter) = comp_obs(i, 2)*180.0/3.14159;
  yy4_comp_obs(counter) = comp_obs(i, 3)*180.0/3.14159;
end

 subplot('position', [0.1,0.2,0.8,0.5])
plot(xx4, yy4, 'b.', xx4_comp_obs, yy4_comp_obs, 'r.');
axis([180 320 0 80])            % display only the region seleceted.


ylabel('Latitude (N)', 'fontsize', 12)
xlabel('Longitude (E)', 'fontsize', 12)
title('GPS RO profile locations in CONUS domain, Jan 1, 2003', 'fontsize', 14)

 subplot('position', [0.0,0.0,0.06,0.06])
 text(0.15, 0.15, 'Fig. 1b')
print -dpsc gpsloc.ps

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

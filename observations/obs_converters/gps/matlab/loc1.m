%% loc1

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

load seleobs.july;
err = seleobs(:, 1);
otype = seleobs(:, 9);

counter = 0;
for i = 1:size(seleobs, 1)
  counter = counter + 1;
  xx4(counter) = seleobs(i, 2)*180.0/3.14159;
  yy4(counter) = seleobs(i, 3)*180.0/3.14159;
end

 subplot('position', [0.1,0.3,0.8,0.5])
plot(xx4, yy4, 'b.');
axis([260 300 0 35])            % display only the region seleceted.
ylabel('Latitude (N)', 'fontsize', 12)
xlabel('Longitude (E)', 'fontsize', 12)
title('GPS locations in CONUS domain, Jan 1, 2003', 'fontsize', 14)

 subplot('position', [0.0,0.0,0.06,0.06])
 text(0.15, 0.15, 'Fig. 1b')
print -dpsc seleobs.ps

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

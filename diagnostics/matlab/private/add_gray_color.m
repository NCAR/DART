function [clim, datarange] = add_gray_color(datarange)
%% Augment the colormap and the CLim so that the lowest color index can be forced to a light gray without compromising the data range.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

bob      = colormap;
ncolors  = length(bob);
bob(1,:) = 0.7; % set lowest color to be light gray.
colormap(bob);

if (nargin < 1) 
   datarange = get(gca,'CLim');
end
cmin = datarange(1);
cmax = datarange(2);

if (cmin ~= cmax)
    dz      = linspace(double(cmin),double(cmax),ncolors-1); % desired dynamic range
    dclim   = dz(2) - dz(1);
    newcmin = double(cmin) - dclim;
    clim    = [newcmin cmax]; % add extra bin at bottom that no data will use.
else
    newcmin = double(cmin);  % newcmin must be same type as 'elev'.
    clim    = [double(newcmin-1.5) double(newcmin+1.5)];
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

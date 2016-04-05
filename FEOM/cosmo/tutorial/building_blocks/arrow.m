function h = arrow(start,stop,scale)
%%  ARROW(start,stop,scale)  draw a line with an arrow pointing from start to stop
%  Draw a line with an  arrow at the end of a line
%  start is the x,y point where the line starts
%  stop is the x,y point where the line stops
%  Scale is an optional argument that will scale the size of the arrows
%  It is assumed that the axis limits are already set

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if nargin==2
  xl = get(gca,'xlim');
  yl = get(gca,'ylim');
  xd = xl(2)-xl(1);        % this sets the scale for the arrow size
  yd = yl(2)-yl(1);        % thus enabling the arrow to appear in correct 
  scale = (xd + yd) / 2;   % proportion to the current axis
end

hold on
axis(axis)

xdif = stop(1) - start(1);
ydif = stop(2) - start(2);

theta = atan(ydif/xdif);  % the angle has to point according to the slope
theta
if(xdif>=0)
  scale = -scale;
end

xx = [start(1), stop(1),(stop(1)+0.02*scale*cos(theta+pi/16)),NaN,stop(1),... 
(stop(1)+0.02*scale*cos(theta-pi/16))]';
yy = [start(2), stop(2), (stop(2)+0.02*scale*sin(theta+pi/16)),NaN,stop(2),... 
(stop(2)+0.02*scale*sin(theta-pi/16))]';

h = plot(xx,yy);
set(h,'LineWidth',2.0);

hold off

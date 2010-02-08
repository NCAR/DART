function h = ray(x,y)
%% DART:ray  plot a single ray/arrow between a pair of points 
% if X is a matrix, each row defines a new set of
% vectors, each column pertains to the beginning and
% ending points for the ray.
%
%
% x = [1 2; 4 3; 5 7; 7 5; 1 1; 3 3];
% y = [9 9; 8 8; 2 4; 7 5; 4 5; 4 2];
% clf; ray(x,y); axis([0 10 0 10]); grid

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (size(x) ~= size(y))
   error('X and Y must be same size') 
end

[m,n] = size(x);

if ( n > 2 )
   error('Can only be two columns [x1 x2] ...') 
end

X = x(:,1);
Y = y(:,1);
U = x(:,2)-x(:,1);
V = y(:,2)-y(:,1);

% use a (private) modified quiver plot ...

quiver(X,Y,U,V,0);


function hh = quiver(varargin)
%QUIVER Quiver plot.
%   QUIVER(X,Y,U,V) plots velocity vectors as arrows with components (u,v)
%   at the points (x,y).  The matrices X,Y,U,V must all be the same size
%   and contain corresponding position and velocity components (X and Y
%   can also be vectors to specify a uniform grid).  QUIVER automatically
%   scales the arrows to fit within the grid.
%
%   QUIVER(U,V) plots velocity vectors at equally spaced points in
%   the x-y plane.
%
%   QUIVER(U,V,S) or QUIVER(X,Y,U,V,S) automatically scales the 
%   arrows to fit within the grid and then stretches them by S.  Use
%   S=0 to plot the arrows without the automatic scaling.
%
%   QUIVER(...,LINESPEC) uses the plot linestyle specified for
%   the velocity vectors.  Any marker in LINESPEC is drawn at the base
%   instead of an arrow on the tip.  Use a marker of '.' to specify
%   no marker at all.  See PLOT for other possibilities.
%
%   QUIVER(...,'filled') fills any markers specified.
%
%   H = QUIVER(...) returns a vector of line handles.
%
%   Example:
%      [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%      z = x .* exp(-x.^2 - y.^2); [px,py] = gradient(z,.2,.15);
%      contour(x,y,z), hold on
%      quiver(x,y,px,py), hold off, axis image
%
%   See also FEATHER, QUIVER3, PLOT.

% Arrow head parameters

alpha      = 0.33; % Size of arrow head relative to the length of the vector
beta       = 0.33; % Width of the base of the arrow head relative to the length
autoscale  = 1;    % Autoscale if ~= 0 then scale by this.
plotarrows = 1;    % Plot arrows
sym = '';

filled = 0;
ls = '-';
ms = '';
col = '';

nin = nargin;
% Parse the string inputs
while isstr(varargin{nin}),
  vv = varargin{nin};
  if ~isempty(vv) & strcmp(lower(vv(1)),'f')
    filled = 1;
    nin = nin-1;
  else
    [l,c,m,msg] = colstyle(vv);
    if ~isempty(msg), 
      error(sprintf('Unknown option "%s".',vv));
    end
    if ~isempty(l), ls = l; end
    if ~isempty(c), col = c; end
    if ~isempty(m), ms = m; plotarrows = 0; end
    if isequal(m,'.'), ms = ''; end % Don't plot '.'
    nin = nin-1;
  end
end

error(nargchk(2,5,nin));

% Check numeric input arguments
if nin<4, % quiver(u,v) or quiver(u,v,s)
  [msg,x,y,u,v] = xyzchk(varargin{1:2});
else
  [msg,x,y,u,v] = xyzchk(varargin{1:4});
end
if ~isempty(msg), error(msg); end

if nin==3 | nin==5, % quiver(u,v,s) or quiver(x,y,u,v,s)
  autoscale = varargin{nin};
end

% Scalar expand u,v
if prod(size(u))==1, u = u(ones(size(x))); end
if prod(size(v))==1, v = v(ones(size(u))); end

if autoscale,
  % Base autoscale value on average spacing in the x and y
  % directions.  Estimate number of points in each direction as
  % either the size of the input arrays or the effective square
  % spacing if x and y are vectors.
  if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
  delx = diff([min(x(:)) max(x(:))])/n;
  dely = diff([min(y(:)) max(y(:))])/m;
  del = delx.^2 + dely.^2;
  if del>0
    len = sqrt((u.^2 + v.^2)/del);
    maxlen = max(len(:));
  else
    maxlen = 0;
  end
  
  if maxlen>0
    autoscale = autoscale*0.9 / maxlen;
  else
    autoscale = autoscale*0.9;
  end
  u = u*autoscale; v = v*autoscale;
end

ax = newplot;
next = lower(get(ax,'NextPlot'));
hold_state = ishold;

% Make velocity vectors
x = x(:).'; y = y(:).';
u = u(:).'; v = v(:).';
uu = [x;x+u;repmat(NaN,size(u))];
vv = [y;y+v;repmat(NaN,size(u))];

h1 = plot(uu(:),vv(:),[col ls]);

if plotarrows,
  % Make arrow heads and plot them
% hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
%       x+u-alpha*(u-beta*(v+eps));repmat(NaN,size(u))];
% hv = [y+v-alpha*(v-beta*(u+eps));y+v; ...
%       y+v-alpha*(v+beta*(u+eps));repmat(NaN,size(v))];

  un = u ./ sqrt(u.^2 + v.^2);   % TJH ... normalize
  vn = v ./ sqrt(u.^2 + v.^2);   % TJH ... normalize

  hu = [x+u- alpha*(un+beta*(vn+eps))  ;x+u; ...
        x+u- alpha*(un-beta*(vn+eps))  ;repmat(NaN,size(u))];
  hv = [y+v- alpha*(vn-beta*(un+eps))  ;y+v; ...
        y+v- alpha*(vn+beta*(un+eps))  ;repmat(NaN,size(v))];

  hold on;
  h2 = plot(hu(:),hv(:),[col ls]);

else
  h2 = [];
end

if ~isempty(ms), % Plot marker on base
  hu = x; hv = y;
  hold on
  h3 = plot(hu(:),hv(:),[col ms]);
  if filled, set(h3,'markerfacecolor',get(h1,'color')); end
else
  h3 = [];
end

if ~hold_state, hold off, view(2); set(ax,'NextPlot',next); end

if nargout>0, hh = [h1;h2;h3]; end

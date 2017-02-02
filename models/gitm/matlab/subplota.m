function ax = subplota(r,c,ri,ci,d,o)
%ax = SUBPLOTA(r,c,ri,ci,d,o)
% r is number of rows
% c is number of columns
% ri is row index
% ci is col index (set to NaN to treat ri as a linear index (go from SW to NE corner of the figure, E first, N second))
% d is distance from edge of "drawing area" to axis (and btw axes)
% o is offset from lower left corner of the figure to the "drawing area"
%
% Example 1: 
% subplota(2,3,2,nan,.01,.01)
%
% Ex2:
% for i=1:6; subplota(2,3,i,nan,.01,.1); plot(1:i,sin(1:i)); end

% DART $Id$
% CREDIT: Alexey Morozov

if isnan(ci)
    i=ri;
    ci=mod(i-1,c)+1;
    ri=floor((i-1)/c)+1;
end

h  = (1-o-(r+1)*d)/r; %height of each axis
w  = (1-o-(c+1)*d)/c; %width of each axis
sc = o+ci*d+(ci-1)*w; %start of each axis in y
sr = o+ri*d+(ri-1)*h; %start of each axis in x

p=[sc sr w h];

ax = axes('Position', p );

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

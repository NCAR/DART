function [face, len] = get_face(lat, lon)
% Routine to efficiently interpolate on aether cube sphere grid
% Returns which face contains (lat, lon) and the length from the edge of the point
% along each of the axes.

% Range adjustment
if(lon == 2*pi) lon = 0; end

% Convert lat lon to x y z on unit sphere
x = cos(lat) * cos(lon);
y = cos(lat) * sin(lon);
z = sin(lat);
vec = [x; y; z];


% Get the longitudes for this point in the two rotated spaces

%====================================================================
% Following code shows individual rotations;
% Single matrix at the end multiplies them together off-line for single rotation
% Can collapse these to a single rotation vector for efficiency
% Rotation 90 degrees around y to put pole on equator
%RY = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
% Rotate 45 degrees around x
%RX = [1 0 0; 0 cosd(45) sind(45); 0 -sind(45) cosd(45)];
% Then rotate 45 degrees around Z
%RZ = [cosd(45) sind(45) 0; -sind(45) cosd(45) 0; 0 0 1];
% Get longitude in the two rotated spaces
%rot_vec = RZ * RX * RY * vec;
%====================================================================
inv_sqrt_2 = 1 / sqrt(2);
rot_vec = [0.5 0.5 -inv_sqrt_2; 0.5 0.5 inv_sqrt_2; inv_sqrt_2 -inv_sqrt_2 0] * vec;
% Compute the longitude in the rotated space
rlon = atan2(rot_vec(2, :), rot_vec(1, :));
if(rlon < 0) rlon = rlon + 2*pi; end

%====================================================================
% Can collapse these to a single rotation vector for efficiency
% Rotation 90 degrees around y to put pole on equator
%RY = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
% Rotate -45 degrees around x
%RX2 = [1 0 0; 0 cosd(-45) sind(-45); 0 -sind(-45) cosd(-45)];
% Then rotate 45 degrees around Z
%RZ = [cosd(45) sind(45) 0; -sind(45) cosd(45) 0; 0 0 1];
% Get longitude in the two rotated spaces
%rot_vec2 = RZ * RX2 * RY * vec;
%====================================================================
rot_vec2 = [-0.5 0.5 -inv_sqrt_2; -0.5 0.5 inv_sqrt_2; inv_sqrt_2 inv_sqrt_2 0] * vec;
% Compute the longitude in the rotated space
rlon2 = atan2(rot_vec2(2, :), rot_vec2(1, :));
if(rlon2 < 0) rlon2 = rlon2 + 2*pi; end



% Which side could we be on, 1 to 4
side = floor(lon / (pi/2)) + 1;
% Which rotated 1 side are we on 
rside = floor(rlon / (pi/2)) + 1;
% Which rotated 2 side
rside2 = floor(rlon2 / (pi/2)) + 1;

% Figure out the face from here (0 to 5, 4 is south, 5 is north)
% These are consistent with the numbering on Aether grid files for the cubed sphere
if(side == 1 && rside == 1)       face = 0; lon_grid(1) = lon;  lon_grid(2) = rlon;
elseif(side == 2 && rside2 == 1)  face = 1; lon_grid(1) = lon;  lon_grid(2) = rlon2;
elseif(side == 3 && rside == 3)   face = 2; lon_grid(1) = lon;  lon_grid(2) = rlon;
elseif(side == 4 && rside2 == 3)  face = 3; lon_grid(1) = lon;  lon_grid(2) = rlon2;
elseif(rside == 4 && rside2 == 4) face = 4; lon_grid(1) = rlon; lon_grid(2) = rlon2;
elseif(rside == 2 && rside2 == 2) face = 5; lon_grid(1) = rlon; lon_grid(2) = rlon2;
else
   [lat, lon, rlon]
   [side, rside, rside2]
   stop
end

% Can also use the fact that the projection is equidistant on the imbedded cube to get what fraction 
% across the imbedded rectangle we are
% Take the longitudes and turn them into a number between -sqrt(1/3) and sqrt(1/3)
lon_grid_m = mod(lon_grid, pi/2);

% Use law of sines to go from lon back to position along edge of imbedded cube 
% The triangle of interest has a side of length 2sqrt(1/3) (1/2 of the planar diagonal of the imbedded cube)
% The angles adjacent to this side are the longitude and 45 degrees
% The angle opposite the side of length 2sqrt(1/3) is pi - (longitude + pi/4)
% The side opposite the longitude is how far along the side of the cube
% The cube side is 2sqrt(1/3), so the length along the side is between zero and this value
gama = pi - (pi/4 + lon_grid_m(1));
len(1) = sqrt(2/3) * sin(lon_grid_m(1)) / sin(gama);

gamb = pi - (pi/4 + lon_grid_m(2));
len(2) = sqrt(2/3) * sin(lon_grid_m(2)) / sin(gamb);

% If we are on sides 2 or 3, the lengths need to be modified because the grid storage
% for Aether goes from smallest latitude to largest and the longitudes of the shifted
% poles are going the opposite way
if(face == 2 | face == 3) len(2) = 2 * sqrt(1/3) - len(2); end

% Same for face 4 (the bottom) but it's the other coordinate that's reversed
if(face == 4) len(1) = 2*sqrt(1/3) - len(1); end








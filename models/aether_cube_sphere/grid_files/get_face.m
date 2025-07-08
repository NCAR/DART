function [face] = get_face(latd, lond)
% Routine to efficiently interpolate on aether cube sphere grid

% Range adjustment
if(lond == 360) lond = 0; end

% Convert lat lon to x y z on unit sphere
x = cosd(latd) * cosd(lond);
y = cosd(latd) * sind(lond);
z = sind(latd);
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
rlond = atan2d(rot_vec(2, :), rot_vec(1, :));
if(rlond < 0) rlond = rlond + 360; end

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
rlond2 = atan2d(rot_vec2(2, :), rot_vec2(1, :));
if(rlond2 < 0) rlond2 = rlond2 + 360; end



% Which side could we be on, 1 to 4
side = floor(lond / 90) + 1;
% Which rotated 1 side are we on 
rside = floor(rlond / 90) + 1;
% Which rotated 2 side
rside2 = floor(rlond2 / 90) + 1;

% Figure out the face from here (1 to 6, 5 is north, 6 is south)
if(side == 1 && rside == 1) face = 1;
elseif(side == 2 && rside2 == 1) face = 2;
elseif(side == 3 && rside == 3) face = 3;
elseif(side == 4 && rside2 == 3) face = 4;
elseif(rside == 2 && rside2 == 2) face = 5;
elseif(rside == 4 && rside2 == 4) face = 6;
else
   [latd, lond, rlond]
   [side, rside, rside2]
   stop
end


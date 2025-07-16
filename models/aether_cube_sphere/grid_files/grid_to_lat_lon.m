% Given the face (from 0 to 5) the number of lons/lats on a face (np) and
% the indices of the i and j grid point, returns the latitude and longitude of the point
function [lat, lon] = grid_to_lat_lon(face, lat_ind, lon_ind, np)

% Some geometry about the grid point distribution
% Cube side is divided into np-1 interior intervals of width 2sqrt(1/3) / np and
% two exterior intervals of half  width, sqrt(1/3) / np 
cube_side = 2 * sqrt(1 / 3);
del = cube_side / np;
half_del = del / 2;
x = sqrt(1/3) - (half_del + del * (lon_ind - 1));
if(face == 5)
   y = sqrt(1/3) - (half_del + del * (lat_ind - 1));
else
   y = -sqrt(1/3) + (half_del + del * (lat_ind - 1));
end

if(face < 4)
   blon = atan2(sqrt(1 / 3), x);
   blat = atan2(y, sqrt(1/3 + x^2));
   blon = blon - pi/4;

   % Above is for face 0; add pi/2 for each additional face tangent to equator
   lon = blon + pi/2 * face; 
   lat = blat;
elseif(face == 4 | face == 5)
   % Face 4 is tangent to south pole
   lon = atan2(y, x);
   lat = atan2(sqrt(1/3), sqrt(x^2 + y^2));

   if(face == 4) lat = -lat; end

   %  Get ready for rotation
   xt = cos(lat) * cos(lon);
   yt = cos(lat) * sin(lon);
   zt = sin(lat);
   vect = [xt; yt; zt];

   % Then rotate 45 degrees around Z
   rot_angle = -pi/4;
   RZ = [cos(rot_angle) sin(rot_angle) 0; -sin(rot_angle) cos(rot_angle) 0; 0 0 1];
   rot_vect = RZ * vect;

   lat = asin(rot_vect(3));
   lon = atan2(rot_vect(2), rot_vect(1));
   % Note that there are inconsistent treatments of the value near longitude
   % 0 in the grid files for Aether. Some points have a value near or just less
   % than 2pi, other points have values just greater than 0. This code tries
   % to avoid values near to 2pi and have 0 instead.
   if(lon < 0) lon =  lon + 2*pi; end
   if(lon >= 2*pi) lon = 0; end

end

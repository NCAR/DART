% Initial work on finding the actual 3 or 4 grid points that create triangle or quad that
% contain a given lat/lon on the Aether cubed sphere.

% Parameter that has the number of model grid points along each dimension 
% This does not include halos; the points are offset from the boundaries
np = 18;

% Some geometry about the grid point distribution
% Cube side is divided into np-1 interior intervals of width 2sqrt(1/3) / np and
% two exterior intervals of half  width, sqrt(1/3) / np 
cube_side = 2 * sqrt(1 / 3); 
del = cube_side / np;
half_del = del / 2;

% Precompute the positions of the points along the cube
% These are only used to verify that these grid computations give the same values as the grid files
for i = 1:18
   x(i) = -sqrt(1/3) + half_del + del * (i - 1);
   box_lon(i) = atan2(sqrt(1 / 3), x(i));
end
box_lon = box_lon - pi / 4;

%-------------------------------------------
% Plot a background spherical surface
th = linspace(0,2*pi, 1000) ;
phi = linspace(-pi/2, pi/2, 1000) ;
[T,P] = meshgrid(th,phi) ;

X = cos(P) .* cos(T);
Y = cos(P) .* sin(T);
Z = sin(P);

% Make the grey scale vary by face
for i = 1:size(X, 1)
   for j = 1:size(X, 2)
      lat(i, j) = asin(Z(i, j));
      lon(i, j) = atan2(Y(i, j), X(i, j));
      if(lon(i, j) < 0) lon(i, j) =  lon(i, j) + 2*pi; end
      if(lon(i, j) == 2*pi) lon(i, j) = 0; end
      [face, lon_grid, len] = get_face(lat(i, j), lon(i, j));
      C(i, j) = face;
   end
end

h = surf(X,Y,Z, C);
set(h, 'linestyle', 'none')
axis equal
% Change the color map
for i = 1:size(colormap, 1)
   c_grey(i, 1:3) = 0.6 + (i / size(colormap, 1) * 0.2);
end
colormap(c_grey);
colorbar
hold on
%-------------------------------------------

% Open grid files for each face for comparison
for face = 0:5
   % Generate the grid file name
   fname = strcat('grid_g000', num2str(face), '.nc');
   xt = ncread(fname, 'Longitude');
   glon(face + 1, :, :) = squeeze(xt(1, 3:end-2, 3:end-2));
   yt = ncread(fname, 'Latitude');
   glat(face + 1, :, :) = squeeze(yt(1, 3:end-2, 3:end-2));

   % Get x, y, z for this face for plotting
   x = cos(glat(face + 1, :, :)) .* cos(glon(face + 1, :, :));
   y = cos(glat(face + 1, :, :)) .* sin(glon(face + 1, :, :));
   z = sin(glat(face + 1, :, :));
   % Plot points on each face with distinct colors 
   h = plot3(squeeze(x), squeeze(y), squeeze(z), '.', 'markersize', 24, 'color', 'blue', 'linewidth', 16);
   hold on

end
%-------------------------------------------

% Enter lats and lons in degrees for starters
pt_lon_d = 1;
pt_lat_d = 0;
% Convert to radians since this is generally used in algorithm
pt_lon = deg2rad(pt_lon_d);
pt_lat = deg2rad(pt_lat_d);
% Convert to x, y, z for plotting
px = cos(pt_lat) .* cos(pt_lon);
py = cos(pt_lat) .* sin(pt_lon);
pz = sin(pt_lat);

% Get the face, the longitudes on the two projections that don't have a pole in this face
% and the length along the two imbedded cube faces.
[face, lon_grid, len] = get_face(pt_lat, pt_lon); 

% Figure out which interval this is in along each cube face; This gives 0 to np intervals
low_grid = floor((len + half_del) / del);
hi_grid = low_grid + 1;

% Plot the point in green
plot3(px, py, pz, '.', 'markersize', 24, 'color', 'green', 'linewidth', 16);

% Straightfoward if not in one of the boundaries
if(all(low_grid > 0) && all(hi_grid < np)) 

   % Highlight the bounding points;
   bp_lon(1) = glon(face + 1, low_grid(2), low_grid(1)); bp_lat(1) = glat(face + 1, low_grid(2), low_grid(1));
   bp_lon(2) = glon(face + 1, low_grid(2), hi_grid(1));  bp_lat(2) = glat(face + 1, low_grid(2), hi_grid(1));
   bp_lon(3) = glon(face + 1, hi_grid(2),  low_grid(1)); bp_lat(3) = glat(face + 1, hi_grid(2),  low_grid(1));
   bp_lon(4) = glon(face + 1, hi_grid(2),  hi_grid(1));  bp_lat(4) = glat(face + 1, hi_grid(2),  hi_grid(1));
else
   % Start with face 1, then generalize
   if(face == 0)
    
      if(low_grid(1) == 0) 
         % On left side
         if(low_grid(2) == 0)
            % Lower left corner
         elseif(hi_grid(2) == np)
            % Upper left Corner
         else
            % Interior left side; Face 3 is to the left
            bp_lon(1) = glon(3 + 1, low_grid(2), np);          bp_lat(1) = glat(3 + 1, low_grid(2), np);
            bp_lon(2) = glon(0 + 1, low_grid(2), hi_grid(1));  bp_lat(2) = glat(0 + 1, low_grid(2), hi_grid(1));
            bp_lon(3) = glon(3 + 1, hi_grid(2),  np);          bp_lat(3) = glat(3 + 1, hi_grid(2), np);
            bp_lon(4) = glon(0 + 1, hi_grid(2),  hi_grid(1));  bp_lat(4) = glat(0 + 1, hi_grid(2),  hi_grid(1));
         end
      elseif(low_grid(1) == np)
         if(low_grid(2) == 0)
            % Lower right corner
         elseif(hi_grid(2) == np)
            % Upper right corner
         else
            % Interior right side; Face 1 is to the right
            bp_lon(1) = glon(0 + 1, low_grid(2), hi_grid(2));  bp_lat(1) = glat(1 + 1, low_grid(2), hi_grid(2));
            bp_lon(2) = glon(1 + 1, low_grid(2), hi_grid(1));  bp_lat(2) = glat(0 + 1, low_grid(2), hi_grid(1));
            bp_lon(3) = glon(1 + 1, hi_grid(2),  np);          bp_lat(3) = glat(1 + 1, hi_grid(2), np);
            bp_lon(4) = glon(0 + 1, hi_grid(2),  hi_grid(1));  bp_lat(4) = glat(0 + 1, hi_grid(2),  np);
         end
      elseif(low_grid(2) == 0)
         % Lower side
      elseif(hi_grid(2) == 0)
         % Top side
      end
   end

end


% Convert each to x, y, z and plot
for i = 1:4
   bpx = cos(bp_lat(i)) .* cos(bp_lon(i));
   bpy = cos(bp_lat(i)) .* sin(bp_lon(i));
   bpz = sin(bp_lat(i));
   % Plot the bounding points in red
   plot3(bpx, bpy, bpz, '.', 'markersize', 24, 'color', 'red', 'linewidth', 16);
end


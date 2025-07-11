r = 1 ;
th = linspace(0,2*pi, 200) ;
phi = linspace(-pi/2, pi/2, 200) ;
[T,P] = meshgrid(th,phi) ;

X = r * cos(P) .* cos(T);
Y = r * cos(P) .* sin(T);
Z = r * sin(P);

% Make the color the face returned by get_face.m
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

surf(X,Y,Z, C);
axis equal
hold on

% Rotation 90 degrees around y to put pole on equator
RY = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
% Rotate 45 degrees around x
RX = [1 0 0; 0 cosd(45) sind(45); 0 -sind(45) cosd(45)];
% Then rotate 45 degrees around Z
RZ = [cosd(45) sind(45) 0; -sind(45) cosd(45) 0; 0 0 1];
% Or rotate z the other way
NRZ = [cosd(135) sind(135) 0; -sind(135) cosd(135) 0; 0 0 1];

% First figure out where to stop each boundary line
% This is where the longitude in the rotated boundary is also pi/2
lon = 0;
lat = -pi/4:0.001:pi/4;
xt = cos(lat) * cos(lon);
yt = cos(lat) * sin(lon);
zt = sin(lat);
vect = [xt; yt; zt];

% Do first rotation
rot_vect = RZ * RX * RY * vect;

% Compute the longitude in the rotated space
rlont = atan2(rot_vect(2, :), rot_vect(1, :));

% Find the index that is closest to zero
[min_val, min_index] = min(abs(rlont));
lat_min = lat(min_index);
[max_val, max_index] = min(abs(rlont - pi/2));
lat_max = lat(max_index);


% Now add on single highlighted longitude great circles
for lon = 0:pi/2:3*pi/2
   % Do a whole bunch of latitudes to get a nice line
   lat = lat_min:0.001:lat_max;
   x = cos(lat) * cos(lon);
   y = cos(lat) * sin(lon);
   z = sin(lat);
   vec = [x; y; z];

   % Do first rotation
   rot_vec = RZ * RX * RY * vec;

   % Do second rotation
   rot_vec2 = NRZ * RX * RY * vec;

   % Plot the three sets of boundary lines over the appropriate range
   plot3(x, y, z, 'color', 'r', 'linewidth', 3);
   plot3(rot_vec(1, :), rot_vec(2, :), rot_vec(3, :), 'color', 'b', 'linewidth', 3);
   plot3(rot_vec2(1, :), rot_vec2(2, :), rot_vec2(3, :), 'color', 'g', 'linewidth', 3);
end

hold on
plot3_grid();

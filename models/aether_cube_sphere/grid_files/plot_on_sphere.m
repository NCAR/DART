r = 1 ;
th = linspace(0,2*pi, 200) ;
phi = linspace(-pi/2, pi/2, 200) ;
[T,P] = meshgrid(th,phi) ;

X = r * cos(P) .* cos(T);
Y = r * cos(P) .* sin(T);
Z = r * sin(P);

% Temporarily make the color the face returned by get_face.m
for i = 1:size(X, 1)
   for j = 1:size(X, 2)
      lat(i, j) = asind(Z(i, j));
      lon(i, j) = atan2d(Y(i, j), X(i, j));
      if(lon(i, j) < 0) lon(i, j) =  lon(i, j) + 360; end
      face = get_face(lat(i, j), lon(i, j));
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

% The first rotation puts the poles at z = 0 and 45 and 225 degrees
% The second rotation puts the poles at z = 0 and 135 and 315 degrees

% First figure out where to stop each boundary line
% This is where the longitude in the rotated boundary is also 90
lond = 0;
latd = -45:0.01:45;
xt = cosd(latd) * cosd(lond);
yt = cosd(latd) * sind(lond);
zt = sind(latd);
vect = [xt; yt; zt];

% Do first rotation
rot_vect = RZ * RX * RY * vect;

% Compute the longitude in the rotated space
rlont = atan2d(rot_vect(2, :), rot_vect(1, :));

% Find the index that is closest to zero
[min_val, min_index] = min(abs(rlont));
lat_min = latd(min_index);
[max_val, max_index] = min(abs(rlont - 90));
lat_max = latd(max_index);


% Now add on single highlighted longitude great circles
for lond = 0:90:270
   % Do a whole bunch of latitudes to get a nice line
   latd = lat_min:0.01:lat_max;
   x = cosd(latd) * cosd(lond);
   y = cosd(latd) * sind(lond);
   z = sind(latd);
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

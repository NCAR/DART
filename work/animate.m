% Matlab routine to do simple read of netcdf file and production of
% animation.  Lots of things could be improved.

% Use netcdf utilities to read in field f
f = getcdf_batch('one_d', 'f')
n = size(f, 2)

% Do plot of first step to set axes
plot(f(:, 1))
axes = axis;
label1 = axes(1) + 0.05 * (axes(2) - axes(1))
label2 = axes(3) + (axes(4) - axes(3)) / 2 

% Initialize movie and loop through each frame
M = moviein(n);
for j = 1:n
   plot(f(:,j))
   axis(axes);
   text(label1, label2, ['Frame ', num2str(j)]);
   M(:,j) = getframe;
end

% Play the movie
movie(M, 5)

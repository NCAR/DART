fid = fopen('long.dat');

a = fscanf(fid,'%f',inf);

field = reshape(a, [89, 89]);

[C,h] = contour ( field([1:20],[80:89])');

clabel(C, h);

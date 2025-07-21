% Uses easy to code 3d method to determine if point is in triangle given
% three d coordinates of coplanar vertices and a point
% If point is inside the triangle, the sum of the angles between the vectors connecting
% the point and the vertices should be 2pi. If the point is outside, the sum is not 2pi.

% This code actually works for vertices in 2 or 3 dimensions. The 3-dimensional case
% has to have coplanar points.

function [inside] = is_point_in_triangle(v1, v2, v3, p)


% Compute the area of the original triangle and the three sub triangles with Heron's formula
% from the 3d lengths. Is sum of the three subs equal to the whole?

% Compute lengths of the sides
len_s1 = norm((v1 - v2));
len_s2 = norm((v3 - v2));
len_s3 = norm((v1 - v3));

% Compute the lengths from the point p
len_p1 = norm((p - v1));
len_p2 = norm((p - v2));
len_p3 = norm((p - v3));

% Area of triangle
at = heron(len_s1, len_s2, len_s3);

% Compute areas of sub triangles
at1 = heron(len_p1, len_p2, len_s1);
at2 = heron(len_p2, len_p3, len_s2);
at3 = heron(len_p3, len_p1, len_s3);

inside = abs(((at1 + at2 + at3) - at) < 1e-8);








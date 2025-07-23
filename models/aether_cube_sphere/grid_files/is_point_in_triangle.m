% Uses easy to code 3d method to determine if point is in triangle given
% three d coordinates of coplanar vertices and a point
% If point is inside the triangle, the sum of the angles between the vectors connecting
% the point and the vertices should be 2pi. If the point is outside, the sum is not 2pi.

% This code actually works for vertices in 2 or 3 dimensions. The 3-dimensional case
% has to have coplanar points.

function [inside, dif_frac ] = is_point_in_triangle(v1, v2, v3, p)

% Compute the area of the original triangle and the three sub triangles with Heron's formula
% from the 3d lengths. Is sum of the three subs equal to the whole?

% Get the projection of the point p onto the plane containing the triangle
% Start by getting perpendicular vector by cross product
perp = cross((v1-v2), v2-v3);
% Get unit vector in direction of perp
unit_perp = perp / norm(perp);
% Projection of vector from v1 to p on the unit perp vector is how much to move to get to plane
offset = dot((p-v1), unit_perp);
p_proj = p - offset * unit_perp;

% Compute lengths of the sides
len_s1 = norm((v1 - v2));
len_s2 = norm((v3 - v2));
len_s3 = norm((v1 - v3));

% Compute the lengths from the point p
len_p1 = norm((p_proj - v1));
len_p2 = norm((p_proj - v2));
len_p3 = norm((p_proj - v3));

% Area of triangle
at = heron(len_s1, len_s2, len_s3);

% Compute areas of sub triangles
at1 = heron(len_p1, len_p2, len_s1);
at2 = heron(len_p2, len_p3, len_s2);
at3 = heron(len_p3, len_p1, len_s3);

area_dif = at1 + at2 + at3 - at;

% Spherical geometry conversion to x,y,z means that the points aren't really co-planar
% in three dimensions so areas don't get very close to zero. Can still approximately check.
% This test has revealed errors in the identification of boxes on the boundaries and the
% corner triangles. The value of 7e-1 prevents failure to allow all tests to proceed
% but will need to be greatly reduced for final testing.
dif_frac = area_dif / at;
inside = abs(dif_frac) < 0.002;

return

